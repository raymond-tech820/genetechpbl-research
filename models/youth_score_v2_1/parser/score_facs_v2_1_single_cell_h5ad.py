#!/usr/bin/env python3
"""Score FACS Youth Score v2.1 from a single-cell raw-count h5ad file.

The frozen v2.1 model is a donor/sample-level pseudobulk model. This wrapper
aggregates raw single-cell counts to sample-level pseudobulk counts and then
applies the frozen scoring formula.

Input requirements:
  - AnnData h5ad with cells x genes raw, non-negative integer-like counts in X,
    raw.X, or a named layer.
  - adata.obs must contain a donor/sample column, e.g. mouse.id or mouse_id.

Output:
  - One row per donor/sample with calibrated Youth Score and coverage metrics.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Iterable
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _read_matrix(adata: ad.AnnData, layer: str | None) -> sp.csr_matrix:
    if layer:
        if layer not in adata.layers:
            raise KeyError(f"Layer {layer!r} not found in h5ad.")
        matrix = adata.layers[layer]
    elif adata.raw is not None:
        matrix = adata.raw.X
    else:
        matrix = adata.X
    matrix = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    if matrix.shape[0] != adata.n_obs:
        raise ValueError("Expected cells x genes matrix; rows must match adata.obs.")
    if matrix.data.size and np.nanmin(matrix.data) < 0:
        raise ValueError("Counts must be raw non-negative values.")
    sample = matrix.data[: min(matrix.data.size, 100_000)]
    if sample.size and not np.allclose(sample, np.round(sample)):
        raise ValueError("Counts do not look integer-like; provide raw counts, not log-normalized data.")
    return matrix


def _gene_names(adata: ad.AnnData, gene_column: str | None) -> list[str]:
    if gene_column:
        if gene_column not in adata.var.columns:
            raise KeyError(f"Gene column {gene_column!r} not found in adata.var.")
        return adata.var[gene_column].astype(str).tolist()
    return adata.var_names.astype(str).tolist()


def pseudobulk_counts(matrix: sp.csr_matrix, sample_ids: Iterable[str]) -> tuple[sp.csr_matrix, list[str], np.ndarray]:
    sample_ids = np.asarray([str(x) for x in sample_ids], dtype=object)
    if len(sample_ids) != matrix.shape[0]:
        raise ValueError("sample_ids length must match number of cells.")
    samples, inverse = np.unique(sample_ids, return_inverse=True)
    design = sp.csr_matrix((np.ones(matrix.shape[0]), (inverse, np.arange(matrix.shape[0]))), shape=(len(samples), matrix.shape[0]))
    bulk = (design @ matrix).tocsr()
    cell_counts = np.asarray(design.sum(axis=1)).ravel().astype(int)
    return bulk, samples.astype(str).tolist(), cell_counts


def score_pseudobulk(bulk_samples_by_genes: sp.csr_matrix, sample_names: list[str], gene_names: list[str], signature_csv: Path, calibration_csv: Path, model: str, pseudocount: float, cell_counts: np.ndarray) -> pd.DataFrame:
    signature_all = pd.read_csv(signature_csv)
    calibration_all = pd.read_csv(calibration_csv)
    signature = signature_all.loc[signature_all["model"] == model].copy()
    calibration = calibration_all.loc[calibration_all["model"] == model].copy()
    if signature.empty:
        raise ValueError(f"No signature rows for model {model!r}.")
    if calibration.empty:
        raise ValueError(f"No calibration row for model {model!r}.")
    gene_to_idx = {gene: i for i, gene in enumerate(gene_names)}
    signature["input_gene_index"] = [gene_to_idx.get(str(gene), -1) for gene in signature["gene"]]
    usable = signature.loc[(signature["input_gene_index"] >= 0) & np.isfinite(signature["training_sd"]) & (signature["training_sd"] > 0)].copy()
    if usable.empty:
        raise ValueError("No usable signature genes found in input h5ad.")
    idx = usable["input_gene_index"].to_numpy(dtype=int)
    counts = bulk_samples_by_genes[:, idx].astype(np.float64).toarray()
    lib = np.asarray(bulk_samples_by_genes.sum(axis=1)).ravel().astype(np.float64)
    if np.any(lib <= 0):
        raise ValueError("At least one pseudobulk sample has zero total counts.")
    expr = np.log2((counts / lib[:, None]) * 1_000_000.0 + pseudocount)
    z = (expr - usable["training_mean"].to_numpy(dtype=float)[None, :]) / usable["training_sd"].to_numpy(dtype=float)[None, :]
    weights = usable["weight"].to_numpy(dtype=float)
    weights[~np.isfinite(weights) | (weights <= 0)] = 1.0
    modules = usable["module"].astype(str).to_numpy()
    def module_score(module: str) -> np.ndarray:
        mask = modules == module
        if not np.any(mask): return np.full(len(sample_names), np.nan)
        return (z[:, mask] @ weights[mask]) / weights[mask].sum()
    young = module_score("young_high")
    old = module_score("old_high")
    raw = young - old
    cal = calibration.iloc[0]
    calibrated = (raw - float(cal["old_reference_center"])) / float(cal["calibration_denominator"])
    weighted_coverage = float(np.abs(usable["weight"]).sum() / np.abs(signature["weight"]).sum())
    return pd.DataFrame({
        "sample_id": sample_names,
        "model": model,
        "n_cells": cell_counts,
        "pseudobulk_library_size": lib,
        "young_module_score": young,
        "old_module_score": old,
        "raw_score": raw,
        "calibrated_score": calibrated,
        "gene_coverage": usable["gene"].nunique() / len(signature),
        "weighted_coverage": weighted_coverage,
    })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--signature-csv", default="models/facs_v2_1_full_data_frozen_signatures_all_models.csv")
    parser.add_argument("--calibration-csv", default="models/facs_v2_1_full_data_frozen_calibration.csv")
    parser.add_argument("--sample-column", default="mouse.id")
    parser.add_argument("--gene-column", default=None)
    parser.add_argument("--layer", default=None)
    parser.add_argument("--model", default="factorial_medium_original")
    parser.add_argument("--pseudocount", type=float, default=1.0)
    args = parser.parse_args()
    adata = ad.read_h5ad(args.input_h5ad)
    if args.sample_column not in adata.obs.columns:
        raise KeyError(f"Sample column {args.sample_column!r} not found in adata.obs.")
    matrix = _read_matrix(adata, args.layer)
    genes = _gene_names(adata, args.gene_column)
    bulk, samples, cell_counts = pseudobulk_counts(matrix, adata.obs[args.sample_column])
    scores = score_pseudobulk(bulk, samples, genes, Path(args.signature_csv), Path(args.calibration_csv), args.model, args.pseudocount, cell_counts)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(output, index=False)
    print(output)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
