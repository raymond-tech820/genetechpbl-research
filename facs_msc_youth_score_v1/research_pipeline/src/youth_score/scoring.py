"""Load cross-validated ensembles and score new sparse expression matrices."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from scipy.special import expit

from .baselines import BaselineFold, _signature_values
from .evaluation import TemperatureScaler
from .features import TokenizerState, encode_named_counts, log_normalize_counts
from .model import GeneTransformer
from .utils import read_json


def align_named_matrix(
    counts: sp.csr_matrix,
    input_gene_names: Sequence[str],
    target_gene_names: Sequence[str],
) -> tuple[sp.csr_matrix, float]:
    input_lookup: dict[str, int] = {}
    for idx, name in enumerate(input_gene_names):
        input_lookup.setdefault(str(name), idx)
    source_indices = []
    target_indices = []
    for target_idx, gene in enumerate(target_gene_names):
        source_idx = input_lookup.get(str(gene))
        if source_idx is not None:
            source_indices.append(source_idx)
            target_indices.append(target_idx)
    if not source_indices:
        return sp.csr_matrix((counts.shape[0], len(target_gene_names)), dtype=np.float32), 0.0
    selected = counts[:, np.asarray(source_indices, dtype=int)].tocoo()
    target_lookup = np.asarray(target_indices, dtype=np.int32)
    aligned = sp.coo_matrix(
        (selected.data, (selected.row, target_lookup[selected.col])),
        shape=(counts.shape[0], len(target_gene_names)),
    ).tocsr()
    return aligned, len(source_indices) / max(len(target_gene_names), 1)


@torch.inference_mode()
def _predict_encoded(
    model: GeneTransformer,
    encoded: dict[str, np.ndarray],
    batch_size: int,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    logits = []
    ages = []
    use_amp = device.type == "cuda"
    for start in range(0, len(encoded["gene_ids"]), batch_size):
        stop = min(start + batch_size, len(encoded["gene_ids"]))
        with torch.amp.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=use_amp):
            batch_logits, batch_age = model(
                torch.from_numpy(encoded["gene_ids"][start:stop].astype(np.int64)).to(device),
                torch.from_numpy(encoded["expression_ids"][start:stop].astype(np.int64)).to(device),
                torch.from_numpy(encoded["rank_ids"][start:stop].astype(np.int64)).to(device),
                torch.from_numpy(encoded["attention_mask"][start:stop]).to(device),
            )
        logits.append(batch_logits.float().cpu().numpy())
        ages.append(batch_age.float().cpu().numpy() * 30.0)
    return np.concatenate(logits), np.concatenate(ages)


@dataclass
class FoldAssets:
    fold: int
    state: TokenizerState
    transformer_models: list[GeneTransformer]
    transformer_temperatures: list[float]
    baseline: BaselineFold
    baseline_temperatures: dict[str, float]


class YouthScoreEnsemble:
    def __init__(self, output_directory: str | Path, device: str = "auto", batch_size: int = 64) -> None:
        self.output_directory = Path(output_directory)
        self.selection = read_json(self.output_directory / "selection.json")
        self.official_model = str(self.selection["official_model"])
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(device)
        self.batch_size = batch_size
        self.folds: list[FoldAssets] = []

        for fold_dir in sorted((self.output_directory / "folds").glob("fold_*")):
            fold = int(fold_dir.name.split("_")[-1])
            state = TokenizerState.load(fold_dir / "tokenizer.json")
            transformer_models = []
            transformer_temperatures = []
            for seed_dir in sorted((fold_dir / "transformer").glob("seed_*")):
                if not (seed_dir / "model.safetensors").exists():
                    continue
                transformer_models.append(GeneTransformer.load(seed_dir, self.device))
                transformer_temperatures.append(float(read_json(seed_dir / "calibration.json")["temperature"]))
            if not transformer_models:
                raise FileNotFoundError(f"No Transformer checkpoints were found in {fold_dir}")
            baseline_dir = fold_dir / "baselines"
            baseline = BaselineFold.load(baseline_dir, state)
            baseline_calibration = read_json(baseline_dir / "baseline_calibration.json")
            baseline_temperatures = {
                key: float(value["temperature"] if isinstance(value, dict) else value)
                for key, value in baseline_calibration.items()
            }
            self.folds.append(
                FoldAssets(
                    fold=fold,
                    state=state,
                    transformer_models=transformer_models,
                    transformer_temperatures=transformer_temperatures,
                    baseline=baseline,
                    baseline_temperatures=baseline_temperatures,
                )
            )
        if not self.folds:
            raise FileNotFoundError(f"No trained folds were found in {self.output_directory}")

    def score_counts(
        self,
        counts: sp.csr_matrix,
        gene_names: Sequence[str],
        cell_ids: Sequence[str] | None = None,
    ) -> pd.DataFrame:
        counts = counts.tocsr()
        if cell_ids is None:
            cell_ids = [str(i) for i in range(counts.shape[0])]
        if len(cell_ids) != counts.shape[0]:
            raise ValueError("The number of cell IDs does not match the expression matrix.")
        official_fold_scores = []
        age_fold_scores = []
        overlaps = []

        for assets in self.folds:
            encoded, overlap = encode_named_counts(counts, gene_names, assets.state)
            overlaps.append(overlap)
            transformer_scores = []
            transformer_ages = []
            for model, temperature in zip(assets.transformer_models, assets.transformer_temperatures):
                logits, ages = _predict_encoded(model, encoded, self.batch_size, self.device)
                transformer_scores.append(TemperatureScaler(temperature).transform(logits))
                transformer_ages.append(ages)
            transformer_score = np.mean(transformer_scores, axis=0)
            age_fold_scores.append(np.mean(transformer_ages, axis=0))

            if self.official_model == "gene_transformer":
                official_fold_scores.append(transformer_score)
                continue

            aligned, baseline_overlap = align_named_matrix(counts, gene_names, assets.state.gene_names)
            overlaps[-1] = min(overlap, baseline_overlap)
            normalized = log_normalize_counts(aligned)
            if self.official_model == "gene_signature":
                values = _signature_values(normalized, assets.baseline.young_local, assets.baseline.old_local)  # type: ignore[arg-type]
                logits = assets.baseline.signature_model.decision_function(values.reshape(-1, 1))  # type: ignore[union-attr]
            elif self.official_model == "elastic_net":
                logits = assets.baseline.elastic_model.decision_function(normalized)  # type: ignore[union-attr]
            else:
                raise ValueError(f"Unsupported official model for external scoring: {self.official_model}")
            temperature = assets.baseline_temperatures[self.official_model]
            official_fold_scores.append(expit(logits / temperature))

        youth_score = np.mean(official_fold_scores, axis=0)
        predicted_age = np.mean(age_fold_scores, axis=0)
        minimum_overlap = float(min(overlaps))
        return pd.DataFrame(
            {
                "cell_id": list(cell_ids),
                "youth_score": youth_score,
                "predicted_age_months": predicted_age,
                "model_id": self.official_model,
                "gene_overlap": minimum_overlap,
                "qc_status": np.where(minimum_overlap >= 0.70, "pass", "insufficient_gene_overlap"),
            }
        )
