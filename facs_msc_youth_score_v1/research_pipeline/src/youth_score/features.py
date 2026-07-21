"""Leakage-safe feature selection and rank-based gene tokenization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import scipy.sparse as sp

from .utils import read_json, write_json


def log_normalize_counts(counts: sp.csr_matrix, target_sum: float = 10_000.0) -> sp.csr_matrix:
    matrix = counts.astype(np.float32, copy=True).tocsr()
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    scales = np.divide(target_sum, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    matrix = sp.diags(scales, format="csr") @ matrix
    np.log1p(matrix.data, out=matrix.data)
    matrix.eliminate_zeros()
    return matrix.tocsr()


def select_highly_variable_genes(
    train_counts: sp.csr_matrix,
    eligible: np.ndarray,
    feature_count: int,
    minimum_detection: float = 0.02,
) -> tuple[np.ndarray, sp.csr_matrix]:
    normalized = log_normalize_counts(train_counts)
    n_cells = normalized.shape[0]
    detected = np.asarray((normalized > 0).sum(axis=0)).ravel() / max(n_cells, 1)
    mean = np.asarray(normalized.mean(axis=0)).ravel()
    squared_mean = np.asarray(normalized.power(2).mean(axis=0)).ravel()
    variance = np.maximum(squared_mean - mean**2, 0.0)
    valid = eligible.astype(bool) & (detected >= minimum_detection) & np.isfinite(variance)
    candidate_indices = np.flatnonzero(valid)
    if len(candidate_indices) < feature_count:
        candidate_indices = np.flatnonzero(eligible.astype(bool) & np.isfinite(variance))
    order = np.argsort(variance[candidate_indices], kind="stable")[::-1]
    selected = candidate_indices[order[: min(feature_count, len(order))]].astype(np.int32)
    if len(selected) < min(256, feature_count):
        raise ValueError(f"Too few eligible genes were found: {len(selected)}")
    return selected, normalized


def fit_expression_bins(
    normalized_train: sp.csr_matrix,
    feature_indices: np.ndarray,
    bin_count: int,
    seed: int,
    sample_size: int = 1_000_000,
) -> np.ndarray:
    values = normalized_train[:, feature_indices].data
    if values.size == 0:
        raise ValueError("Cannot fit expression bins without non-zero values.")
    if values.size > sample_size:
        rng = np.random.default_rng(seed)
        values = values[rng.choice(values.size, size=sample_size, replace=False)]
    quantiles = np.linspace(0.0, 1.0, bin_count + 1)[1:-1]
    edges = np.unique(np.quantile(values, quantiles)).astype(np.float32)
    return edges


@dataclass(frozen=True)
class TokenizerState:
    gene_names: tuple[str, ...]
    source_gene_indices: tuple[int, ...]
    expression_edges: tuple[float, ...]
    sequence_length: int
    target_sum: float = 10_000.0

    @property
    def pad_token_id(self) -> int:
        return 0

    @property
    def cls_token_id(self) -> int:
        return len(self.gene_names) + 1

    @property
    def vocabulary_size(self) -> int:
        return len(self.gene_names) + 2

    @property
    def expression_bin_count(self) -> int:
        return len(self.expression_edges) + 1

    def to_dict(self) -> dict:
        return {
            "gene_names": list(self.gene_names),
            "source_gene_indices": list(self.source_gene_indices),
            "expression_edges": list(self.expression_edges),
            "sequence_length": self.sequence_length,
            "target_sum": self.target_sum,
        }

    @classmethod
    def from_dict(cls, value: dict) -> "TokenizerState":
        return cls(
            gene_names=tuple(value["gene_names"]),
            source_gene_indices=tuple(int(v) for v in value["source_gene_indices"]),
            expression_edges=tuple(float(v) for v in value["expression_edges"]),
            sequence_length=int(value["sequence_length"]),
            target_sum=float(value.get("target_sum", 10_000.0)),
        )

    def save(self, path: str | Path) -> None:
        write_json(path, self.to_dict())

    @classmethod
    def load(cls, path: str | Path) -> "TokenizerState":
        return cls.from_dict(read_json(path))


def build_tokenizer_state(
    train_counts: sp.csr_matrix,
    all_gene_names: Sequence[str],
    eligible: np.ndarray,
    feature_count: int,
    sequence_length: int,
    expression_bins: int,
    seed: int,
) -> TokenizerState:
    feature_indices, normalized = select_highly_variable_genes(
        train_counts=train_counts,
        eligible=eligible,
        feature_count=feature_count,
    )
    edges = fit_expression_bins(normalized, feature_indices, expression_bins, seed)
    return TokenizerState(
        gene_names=tuple(str(all_gene_names[i]) for i in feature_indices),
        source_gene_indices=tuple(int(i) for i in feature_indices),
        expression_edges=tuple(float(v) for v in edges),
        sequence_length=sequence_length,
    )


def _encode_rows(
    counts: sp.csr_matrix,
    source_to_token: np.ndarray,
    state: TokenizerState,
) -> dict[str, np.ndarray]:
    counts = counts.tocsr()
    n_cells = counts.shape[0]
    width = state.sequence_length + 1
    gene_ids = np.zeros((n_cells, width), dtype=np.int32)
    expression_ids = np.zeros((n_cells, width), dtype=np.int16)
    rank_ids = np.zeros((n_cells, width), dtype=np.int16)
    attention_mask = np.zeros((n_cells, width), dtype=bool)
    gene_ids[:, 0] = state.cls_token_id
    attention_mask[:, 0] = True
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    edges = np.asarray(state.expression_edges, dtype=np.float32)

    for row in range(n_cells):
        start, end = counts.indptr[row], counts.indptr[row + 1]
        source_indices = counts.indices[start:end]
        values = counts.data[start:end].astype(np.float32, copy=False)
        token_ids = source_to_token[source_indices]
        keep = token_ids > 0
        if not keep.any() or totals[row] <= 0:
            continue
        token_ids = token_ids[keep]
        normalized = np.log1p(values[keep] * (state.target_sum / totals[row]))
        if len(normalized) > state.sequence_length:
            top = np.argpartition(normalized, -state.sequence_length)[-state.sequence_length :]
            token_ids = token_ids[top]
            normalized = normalized[top]
        order = np.argsort(normalized, kind="stable")[::-1]
        token_ids = token_ids[order]
        normalized = normalized[order]
        length = len(token_ids)
        gene_ids[row, 1 : length + 1] = token_ids
        expression_ids[row, 1 : length + 1] = np.searchsorted(edges, normalized, side="right") + 1
        rank_ids[row, 1 : length + 1] = np.arange(1, length + 1, dtype=np.int16)
        attention_mask[row, : length + 1] = True

    return {
        "gene_ids": gene_ids,
        "expression_ids": expression_ids,
        "rank_ids": rank_ids,
        "attention_mask": attention_mask,
    }


def encode_tms_counts(counts: sp.csr_matrix, state: TokenizerState) -> dict[str, np.ndarray]:
    source_to_token = np.zeros(counts.shape[1], dtype=np.int32)
    source_indices = np.asarray(state.source_gene_indices, dtype=np.int64)
    if source_indices.max(initial=-1) >= counts.shape[1]:
        raise ValueError("Tokenizer source gene indices do not fit the input matrix.")
    source_to_token[source_indices] = np.arange(1, len(source_indices) + 1, dtype=np.int32)
    return _encode_rows(counts, source_to_token, state)


def encode_named_counts(
    counts: sp.csr_matrix,
    input_gene_names: Sequence[str],
    state: TokenizerState,
) -> tuple[dict[str, np.ndarray], float]:
    token_lookup = {gene: i + 1 for i, gene in enumerate(state.gene_names)}
    source_to_token = np.zeros(len(input_gene_names), dtype=np.int32)
    matched = 0
    for source_idx, gene in enumerate(input_gene_names):
        token_id = token_lookup.get(str(gene), 0)
        source_to_token[source_idx] = token_id
        matched += int(token_id > 0)
    overlap = matched / max(len(state.gene_names), 1)
    return _encode_rows(counts, source_to_token, state), overlap


def save_encoded(path: str | Path, encoded: dict[str, np.ndarray]) -> None:
    np.savez_compressed(path, **encoded)


def load_encoded(path: str | Path) -> dict[str, np.ndarray]:
    with np.load(path) as loaded:
        return {key: loaded[key] for key in loaded.files}

