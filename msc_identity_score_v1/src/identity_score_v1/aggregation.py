"""Correct animal-level aggregation for Identity Score v1."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .scoring import PRIMARY_COLUMN, RANK_COLUMN


def _p95(values: pd.Series) -> float:
    return float(values.quantile(0.95))


def _summary(frame: pd.DataFrame) -> pd.Series:
    return pd.Series(
        {
            "n_cells": int(len(frame)),
            "identity_primary_median": float(frame[PRIMARY_COLUMN].median()),
            "identity_primary_mean": float(frame[PRIMARY_COLUMN].mean()),
            "identity_primary_p95": _p95(frame[PRIMARY_COLUMN]),
            "identity_rank_median": float(frame[RANK_COLUMN].median()),
            "identity_rank_mean": float(frame[RANK_COLUMN].mean()),
            "identity_rank_p95": _p95(frame[RANK_COLUMN]),
        }
    )


def aggregate_identity_scores(
    cells: pd.DataFrame,
    *,
    dataset_id: str,
    species: str,
    cell_population: str,
    unit_prefix: str,
    gene_coverage: float,
    expected_biological_units: int | None = None,
) -> dict[str, pd.DataFrame]:
    required = {
        "cell_id",
        "age_raw",
        "age_group",
        "treatment_raw",
        "analysis_role",
        "animal",
        "known_animal",
        "state",
        "sample",
        "batch",
        PRIMARY_COLUMN,
        RANK_COLUMN,
    }
    missing = required.difference(cells.columns)
    if missing:
        raise ValueError(f"Cell-score table is missing columns: {sorted(missing)}")
    if cells["cell_id"].duplicated().any():
        raise ValueError("Cell IDs must be unique.")
    if cells[[PRIMARY_COLUMN, RANK_COLUMN]].isna().any().any():
        raise ValueError("Identity scores contain missing values.")

    minimal = cells[
        [
            "cell_id",
            "age_raw",
            "age_group",
            "treatment_raw",
            "analysis_role",
            "animal",
            "known_animal",
            "state",
            "sample",
            "batch",
            PRIMARY_COLUMN,
            RANK_COLUMN,
        ]
    ].copy()
    minimal["biological_unit_id"] = np.where(
        minimal["known_animal"],
        unit_prefix
        + "|"
        + minimal["age_group"].astype(str)
        + "|"
        + minimal["treatment_raw"].astype(str)
        + "|animal_"
        + minimal["animal"].astype(str),
        unit_prefix
        + "|"
        + minimal["age_group"].astype(str)
        + "|"
        + minimal["treatment_raw"].astype(str)
        + "|unknown_pool",
    )

    known = minimal.loc[minimal["known_animal"]].copy()
    donor = (
        known.groupby(
            [
                "age_raw",
                "age_group",
                "treatment_raw",
                "analysis_role",
                "animal",
            ],
            observed=True,
            sort=False,
        )
        .apply(_summary, include_groups=False)
        .reset_index()
    )
    donor.insert(
        0,
        "biological_unit_id",
        unit_prefix
        + "|"
        + donor["age_group"].astype(str)
        + "|"
        + donor["treatment_raw"].astype(str)
        + "|animal_"
        + donor["animal"].astype(str),
    )
    donor.insert(1, "dataset_id", dataset_id)
    donor.insert(2, "species", species)
    donor.insert(3, "cell_population", cell_population)
    donor.rename(
        columns={
            "treatment_raw": "exact_treatment_arm",
            "animal": "animal_label",
        },
        inplace=True,
    )
    donor["animal_label_scope"] = "nested_within_age_and_exact_treatment_arm"
    donor["cross_arm_pairing_status"] = "unpaired"
    donor["identity_gene_coverage"] = float(gene_coverage)
    donor["qc_status"] = "pass"
    donor["n_cells"] = donor["n_cells"].astype(int)
    donor = donor.sort_values(
        ["age_group", "exact_treatment_arm", "animal_label"], kind="stable"
    ).reset_index(drop=True)

    if donor["biological_unit_id"].duplicated().any():
        raise AssertionError("Biological-unit IDs are not unique.")
    if int(donor["n_cells"].sum()) != len(known):
        raise AssertionError("Known-cell count changed during aggregation.")
    if expected_biological_units is not None and len(donor) != expected_biological_units:
        raise AssertionError(
            f"Expected {expected_biological_units} biological units, observed {len(donor)}."
        )

    condition = (
        donor.groupby(
            [
                "dataset_id",
                "species",
                "cell_population",
                "age_raw",
                "age_group",
                "exact_treatment_arm",
                "analysis_role",
            ],
            observed=True,
        )
        .agg(
            n_biological_units=("biological_unit_id", "nunique"),
            total_cells=("n_cells", "sum"),
            identity_primary_median_of_donor_medians=(
                "identity_primary_median",
                "median",
            ),
            identity_primary_q25_of_donor_medians=(
                "identity_primary_median",
                lambda x: x.quantile(0.25),
            ),
            identity_primary_q75_of_donor_medians=(
                "identity_primary_median",
                lambda x: x.quantile(0.75),
            ),
            identity_rank_median_of_donor_medians=(
                "identity_rank_median",
                "median",
            ),
            identity_rank_q25_of_donor_medians=(
                "identity_rank_median",
                lambda x: x.quantile(0.25),
            ),
            identity_rank_q75_of_donor_medians=(
                "identity_rank_median",
                lambda x: x.quantile(0.75),
            ),
        )
        .reset_index()
    )
    condition["statistical_note"] = "descriptive; treatment arms are unpaired"
    condition["n_biological_units"] = condition["n_biological_units"].astype(int)
    condition["total_cells"] = condition["total_cells"].astype(int)

    unknown = minimal.loc[~minimal["known_animal"]].copy()
    if unknown.empty:
        unknown_summary = pd.DataFrame()
    else:
        unknown_summary = (
            unknown.groupby(
                ["age_raw", "age_group", "treatment_raw", "analysis_role"],
                observed=True,
                sort=False,
            )
            .apply(_summary, include_groups=False)
            .reset_index()
            .rename(columns={"treatment_raw": "exact_treatment_arm"})
        )
        unknown_summary.insert(0, "dataset_id", dataset_id)
        unknown_summary["n_cells"] = unknown_summary["n_cells"].astype(int)
        unknown_summary["scope"] = "unknown_animal_cells_descriptive_only"
        unknown_summary["inferential_status"] = "excluded_from_donor_inference"

    return {
        "cells": minimal,
        "donors": donor,
        "conditions": condition,
        "unknown": unknown_summary,
    }


def write_aggregated_outputs(
    tables: dict[str, pd.DataFrame], output_directory: str | Path
) -> None:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    tables["cells"].to_parquet(
        output / "identity_cell_scores_minimal.parquet", index=False
    )
    tables["donors"].to_csv(
        output / "identity_donor_by_condition.csv", index=False
    )
    tables["conditions"].to_csv(
        output / "identity_condition_summary.csv", index=False
    )
    tables["unknown"].to_csv(
        output / "identity_unknown_cell_summary.csv", index=False
    )
