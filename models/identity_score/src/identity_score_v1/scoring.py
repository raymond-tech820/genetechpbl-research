"""Exact cell-level scoring algorithm used by Identity Score v1 (2026-08-06)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .resources import IDENTITY_MODULE


EXCLUDED_PREFIXES = ("MT-", "MT.", "RPL", "RPS", "HBA", "HBB")
PRIMARY_COLUMN = f"{IDENTITY_MODULE}_primary"
RANK_COLUMN = f"{IDENTITY_MODULE}_rank"


@dataclass(frozen=True)
class IdentityPlan:
    positive: np.ndarray
    negative: np.ndarray
    primary_weight: np.ndarray
    coverage: float
    expected_genes: int
    observed_genes: int
    status: str
    duplicate_input_symbols: int


def log_normalize_counts(
    counts: sp.spmatrix,
    target_sum: float = 10_000.0,
) -> sp.csr_matrix:
    """Normalize each cell to ``target_sum`` counts and apply log1p."""
    counts = counts.tocsr().astype(np.float64)
    if counts.data.size and np.nanmin(counts.data) < 0:
        raise ValueError("Identity Score expects non-negative raw counts.")
    total = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.divide(target_sum, total, out=np.zeros_like(total), where=total > 0)
    normalized = sp.diags(scale).dot(counts).tocsr()
    normalized.data = np.log1p(normalized.data)
    return normalized.astype(np.float32)


def canonical_gene_lookup(gene_names: Sequence[str]) -> tuple[dict[str, int], int]:
    lookup: dict[str, int] = {}
    duplicates = 0
    for index, symbol in enumerate(gene_names):
        canonical = str(symbol).strip().upper()
        if canonical in lookup:
            duplicates += 1
        else:
            lookup[canonical] = index
    return lookup, duplicates


def _noise_gene(symbol: str) -> bool:
    return str(symbol).strip().upper().startswith(EXCLUDED_PREFIXES)


def _expression_bins(mean_expression: np.ndarray, bins: int) -> np.ndarray:
    order = np.argsort(mean_expression, kind="stable")
    labels = np.empty(len(order), dtype=int)
    for label, group in enumerate(np.array_split(order, min(bins, len(order)))):
        labels[group] = label
    return labels


def prepare_identity_plan(
    gene_names: Sequence[str],
    identity_gene_table: pd.DataFrame,
    all_module_gene_table: pd.DataFrame,
    reference_means: np.ndarray,
    minimum_gene_coverage: float = 0.70,
    n_control_sets: int = 100,
    control_bins: int = 24,
    seed: int = 20260729,
) -> tuple[IdentityPlan, pd.DataFrame]:
    """Create the deterministic primary-score weight vector.

    The background universe excludes every included gene from the complete
    frozen reprogramming-module table. This detail is required to reproduce the
    2026-08-06 GSE176206 scores; excluding only the 11 identity genes would
    produce a different matched background.
    """
    reference_means = np.asarray(reference_means, dtype=float)
    if len(reference_means) != len(gene_names):
        raise ValueError("reference_means must contain one value per input gene.")
    if n_control_sets < 1 or control_bins < 1:
        raise ValueError("n_control_sets and control_bins must both be positive.")

    lookup, duplicate_symbols = canonical_gene_lookup(gene_names)
    members = identity_gene_table.loc[
        identity_gene_table["included"]
        & identity_gene_table["module"].eq(IDENTITY_MODULE)
    ].copy()
    if members.empty:
        raise ValueError("No included msc_identity_core genes were found.")

    mapped = members["gene_symbol"].astype(str).str.upper().map(lookup)
    observed = members.loc[mapped.notna()].copy()
    indices = mapped.loc[mapped.notna()].astype(int).to_numpy()
    directions = observed["direction"].astype(int).to_numpy()
    positive = indices[directions == 1]
    negative = indices[directions == -1]
    coverage = len(observed) / len(members)
    status = (
        "pass"
        if coverage >= minimum_gene_coverage
        else "insufficient_gene_coverage"
    )

    all_signature = {
        str(symbol).strip().upper()
        for symbol in all_module_gene_table.loc[
            all_module_gene_table["included"], "gene_symbol"
        ]
    }
    universe = np.asarray(
        [
            index
            for index, symbol in enumerate(gene_names)
            if str(symbol).strip().upper() not in all_signature
            and not _noise_gene(str(symbol))
        ],
        dtype=int,
    )
    if len(universe) < 100:
        raise ValueError("Fewer than 100 eligible background genes remain.")

    weight = np.zeros(len(gene_names), dtype=np.float32)
    if status == "pass":
        if len(positive):
            weight[positive] += 1.0 / len(positive)
        if len(negative):
            weight[negative] -= 1.0 / len(negative)

        labels = _expression_bins(reference_means, control_bins)
        universe_by_bin = {
            label: universe[labels[universe] == label]
            for label in np.unique(labels)
        }
        axis_indices = set(map(int, indices))
        control_weight = np.zeros(len(gene_names), dtype=np.float64)
        module_rng = np.random.default_rng(seed)
        for _ in range(n_control_sets):
            for gene_index, direction in zip(indices, directions, strict=True):
                candidates = universe_by_bin.get(
                    labels[gene_index], np.asarray([], dtype=int)
                )
                candidates = np.asarray(
                    [
                        candidate
                        for candidate in candidates
                        if int(candidate) not in axis_indices
                    ],
                    dtype=int,
                )
                if len(candidates) == 0:
                    candidates = universe
                selected = int(module_rng.choice(candidates))
                divisor = len(positive) if direction == 1 else len(negative)
                control_weight[selected] += direction / max(divisor, 1)
        weight -= (control_weight / n_control_sets).astype(np.float32)

    plan = IdentityPlan(
        positive=positive,
        negative=negative,
        primary_weight=weight,
        coverage=float(coverage),
        expected_genes=int(len(members)),
        observed_genes=int(len(observed)),
        status=status,
        duplicate_input_symbols=duplicate_symbols,
    )
    qc = pd.DataFrame(
        [
            {
                "module": IDENTITY_MODULE,
                "expected_genes": len(members),
                "observed_genes": len(observed),
                "gene_coverage": coverage,
                "status": status,
                "positive_genes": len(positive),
                "negative_genes": len(negative),
                "control_sets": n_control_sets,
                "control_bins": control_bins,
                "seed": seed,
                "duplicate_input_symbols": duplicate_symbols,
            }
        ]
    )
    return plan, qc


def _rank_score(counts: sp.csr_matrix, plan: IdentityPlan) -> np.ndarray:
    lookup = np.zeros(counts.shape[1], dtype=np.float32)
    lookup[plan.positive] = 1.0
    lookup[plan.negative] = -1.0
    denominator = max(len(plan.positive) + len(plan.negative), 1)
    output = np.zeros(counts.shape[0], dtype=np.float32)
    counts = counts.tocsr()
    for row in range(counts.shape[0]):
        start, stop = counts.indptr[row], counts.indptr[row + 1]
        row_indices = counts.indices[start:stop]
        row_values = counts.data[start:stop]
        if len(row_values) == 0:
            continue
        sorted_values = np.sort(row_values)
        percentiles = (
            np.searchsorted(sorted_values, row_values, side="right")
            / len(sorted_values)
        )
        output[row] = float(percentiles @ lookup[row_indices]) / denominator
    return output


def score_chunk(
    counts: sp.spmatrix,
    plan: IdentityPlan,
    cell_ids: Sequence[str],
    target_sum: float = 10_000.0,
) -> pd.DataFrame:
    if counts.shape[0] != len(cell_ids):
        raise ValueError("One cell_id is required for every count-matrix row.")
    output = pd.DataFrame({"cell_id": list(map(str, cell_ids))})
    if plan.status != "pass":
        output[PRIMARY_COLUMN] = np.nan
        output[RANK_COLUMN] = np.nan
        return output
    counts = counts.tocsr()
    normalized = log_normalize_counts(counts, target_sum=target_sum)
    output[PRIMARY_COLUMN] = np.asarray(
        normalized @ plan.primary_weight
    ).ravel()
    output[RANK_COLUMN] = _rank_score(counts, plan)
    return output
