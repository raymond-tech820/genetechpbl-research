"""Compact gene-token Transformer used by the Youth Score pipeline."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import torch
from safetensors.torch import load_file, save_file
from torch import nn
from torch.utils.data import Dataset

from .utils import read_json, write_json


@dataclass(frozen=True)
class TransformerShape:
    vocabulary_size: int
    expression_bin_count: int
    sequence_length: int
    d_model: int = 128
    heads: int = 8
    layers: int = 2
    feedforward: int = 512
    dropout: float = 0.2
    token_dropout: float = 0.15


class EncodedCellDataset(Dataset):
    def __init__(
        self,
        encoded: dict[str, np.ndarray],
        labels: np.ndarray,
        ages: np.ndarray,
        weights: np.ndarray,
        indices: np.ndarray,
    ) -> None:
        self.encoded = encoded
        self.labels = labels.astype(np.float32)
        self.ages = ages.astype(np.float32)
        self.weights = weights.astype(np.float32)
        self.indices = np.asarray(indices, dtype=np.int64)

    def __len__(self) -> int:
        return len(self.indices)

    def __getitem__(self, item: int) -> dict[str, torch.Tensor]:
        idx = self.indices[item]
        return {
            "gene_ids": torch.from_numpy(self.encoded["gene_ids"][idx].astype(np.int64, copy=False)),
            "expression_ids": torch.from_numpy(self.encoded["expression_ids"][idx].astype(np.int64, copy=False)),
            "rank_ids": torch.from_numpy(self.encoded["rank_ids"][idx].astype(np.int64, copy=False)),
            "attention_mask": torch.from_numpy(self.encoded["attention_mask"][idx]),
            "label": torch.tensor(self.labels[idx], dtype=torch.float32),
            "age": torch.tensor(self.ages[idx], dtype=torch.float32),
            "weight": torch.tensor(self.weights[idx], dtype=torch.float32),
            "index": torch.tensor(idx, dtype=torch.int64),
        }


class GeneTransformer(nn.Module):
    def __init__(self, shape: TransformerShape) -> None:
        super().__init__()
        self.shape = shape
        self.gene_embedding = nn.Embedding(shape.vocabulary_size, shape.d_model, padding_idx=0)
        self.expression_embedding = nn.Embedding(shape.expression_bin_count + 1, shape.d_model, padding_idx=0)
        self.rank_embedding = nn.Embedding(shape.sequence_length + 1, shape.d_model, padding_idx=0)
        self.input_norm = nn.LayerNorm(shape.d_model)
        self.input_dropout = nn.Dropout(shape.dropout)
        layer = nn.TransformerEncoderLayer(
            d_model=shape.d_model,
            nhead=shape.heads,
            dim_feedforward=shape.feedforward,
            dropout=shape.dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(
            layer,
            num_layers=shape.layers,
            norm=nn.LayerNorm(shape.d_model),
            enable_nested_tensor=False,
        )
        self.youth_head = nn.Sequential(nn.LayerNorm(shape.d_model), nn.Linear(shape.d_model, 1))
        self.age_head = nn.Sequential(
            nn.LayerNorm(shape.d_model),
            nn.Linear(shape.d_model, shape.d_model // 2),
            nn.GELU(),
            nn.Dropout(shape.dropout),
            nn.Linear(shape.d_model // 2, 1),
        )

    def _apply_token_dropout(
        self,
        gene_ids: torch.Tensor,
        expression_ids: torch.Tensor,
        rank_ids: torch.Tensor,
        attention_mask: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        if not self.training or self.shape.token_dropout <= 0:
            return gene_ids, expression_ids, rank_ids, attention_mask
        drop = (torch.rand_like(gene_ids, dtype=torch.float32) < self.shape.token_dropout) & attention_mask
        drop[:, 0] = False
        gene_ids = gene_ids.masked_fill(drop, 0)
        expression_ids = expression_ids.masked_fill(drop, 0)
        rank_ids = rank_ids.masked_fill(drop, 0)
        attention_mask = attention_mask & ~drop
        return gene_ids, expression_ids, rank_ids, attention_mask

    def forward(
        self,
        gene_ids: torch.Tensor,
        expression_ids: torch.Tensor,
        rank_ids: torch.Tensor,
        attention_mask: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        gene_ids, expression_ids, rank_ids, attention_mask = self._apply_token_dropout(
            gene_ids, expression_ids, rank_ids, attention_mask
        )
        hidden = (
            self.gene_embedding(gene_ids)
            + self.expression_embedding(expression_ids)
            + self.rank_embedding(rank_ids)
        )
        hidden = self.input_dropout(self.input_norm(hidden))
        hidden = self.encoder(hidden, src_key_padding_mask=~attention_mask)
        cls = hidden[:, 0]
        return self.youth_head(cls).squeeze(-1), self.age_head(cls).squeeze(-1)

    def save(self, directory: str | Path) -> None:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        save_file(self.state_dict(), str(directory / "model.safetensors"))
        write_json(directory / "model_shape.json", asdict(self.shape))

    @classmethod
    def load(cls, directory: str | Path, device: str | torch.device = "cpu") -> "GeneTransformer":
        directory = Path(directory)
        shape = TransformerShape(**read_json(directory / "model_shape.json"))
        model = cls(shape)
        model.load_state_dict(load_file(str(directory / "model.safetensors"), device=str(device)))
        return model.to(device)
