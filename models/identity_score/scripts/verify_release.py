"""Compare a candidate run with the frozen 2026-08-06 reference outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


TABLES = (
    "identity_donor_by_condition.csv",
    "identity_condition_summary.csv",
    "identity_gene_set_qc.csv",
)
SCORE_COLUMNS = ("msc_identity_core_primary", "msc_identity_core_rank")


def compare_table(candidate_path: Path, reference_path: Path) -> dict:
    candidate = pd.read_csv(candidate_path)
    reference = pd.read_csv(reference_path)
    if list(candidate.columns) != list(reference.columns):
        raise AssertionError(f"Column mismatch for {candidate_path.name}.")
    if candidate.shape != reference.shape:
        raise AssertionError(f"Shape mismatch for {candidate_path.name}.")
    numeric = list(reference.select_dtypes(include=["number"]).columns)
    text = [column for column in reference.columns if column not in numeric]
    if text and not candidate[text].fillna("").equals(reference[text].fillna("")):
        raise AssertionError(f"Text-field mismatch for {candidate_path.name}.")
    maximum = 0.0
    if numeric:
        candidate_values = candidate[numeric].to_numpy(dtype=float)
        reference_values = reference[numeric].to_numpy(dtype=float)
        if not np.allclose(candidate_values, reference_values, rtol=0, atol=1e-6, equal_nan=True):
            raise AssertionError(f"Numeric mismatch for {candidate_path.name}.")
        finite = np.isfinite(candidate_values) & np.isfinite(reference_values)
        if finite.any():
            maximum = float(np.max(np.abs(candidate_values[finite] - reference_values[finite])))
    return {"rows": len(candidate), "max_abs_numeric_difference": maximum}


def compare_cell_scores(candidate_path: Path, reference_path: Path) -> dict:
    candidate = pd.read_parquet(candidate_path)[["cell_id", *SCORE_COLUMNS]]
    reference = pd.read_parquet(reference_path)[["cell_id", *SCORE_COLUMNS]]
    merged = candidate.merge(reference, on="cell_id", how="outer", suffixes=("_candidate", "_reference"), indicator=True, validate="one_to_one")
    if not merged["_merge"].eq("both").all():
        raise AssertionError("Candidate and reference cell IDs differ.")
    differences = [
        np.abs(merged[f"{column}_candidate"].to_numpy(dtype=float) - merged[f"{column}_reference"].to_numpy(dtype=float))
        for column in SCORE_COLUMNS
    ]
    maximum = float(np.max(np.column_stack(differences)))
    if maximum > 1e-6:
        raise AssertionError(f"Cell-level score difference {maximum:.9g} exceeds tolerance.")
    return {"cells": len(merged), "max_abs_score_difference": maximum}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-directory", required=True)
    parser.add_argument("--reference-directory", default=str(Path(__file__).resolve().parents[1] / "reference_results"))
    parser.add_argument("--reference-cell-scores")
    args = parser.parse_args()
    candidate = Path(args.candidate_directory).resolve()
    reference = Path(args.reference_directory).resolve()
    result = {name: compare_table(candidate / name, reference / name) for name in TABLES}
    if args.reference_cell_scores:
        result["cell_scores"] = compare_cell_scores(candidate / "identity_cell_scores_minimal.parquet", Path(args.reference_cell_scores).resolve())
    result["status"] = "pass"
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
