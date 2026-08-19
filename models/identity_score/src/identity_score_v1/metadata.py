"""Dataset metadata normalization with explicit treatment-arm preservation."""

from __future__ import annotations

import numpy as np
import pandas as pd


UNKNOWN_LABELS = {"", "unknown", "nan", "none", "unassigned"}


def is_known_animal(values: pd.Series) -> pd.Series:
    return ~values.astype(str).str.strip().str.lower().isin(UNKNOWN_LABELS)


def normalize_metadata(obs: pd.DataFrame, config: dict) -> pd.DataFrame:
    age_column = config.get("age_column", "age")
    treatment_column = config.get("treatment_column", "treatment")
    animal_column = config.get("animal_column", "animal")
    state_column = config.get("state_column", "state")
    sample_column = config.get("sample_column", "sample")
    batch_column = config.get("batch_column", "batch")
    required = {age_column, treatment_column, animal_column, state_column}
    missing = required.difference(obs.columns)
    if missing:
        raise ValueError(f"Input obs is missing columns: {sorted(missing)}")

    output = pd.DataFrame(index=obs.index)
    output["cell_id"] = obs.index.astype(str)
    output["age_raw"] = obs[age_column].astype(str)
    age = output["age_raw"].str.strip().str.lower()
    young_labels = {str(x).strip().lower() for x in config["young_labels"]}
    aged_labels = {str(x).strip().lower() for x in config["aged_labels"]}
    output["age_group"] = np.select(
        [age.isin(young_labels), age.isin(aged_labels)],
        ["young", "aged"],
        default="unmapped",
    )

    output["treatment_raw"] = obs[treatment_column].astype(str)
    sokm_labels = {str(x) for x in config["sokm_labels"]}
    control_labels = {str(x) for x in config["control_labels"]}
    output["analysis_role"] = np.select(
        [
            output["treatment_raw"].isin(sokm_labels),
            output["treatment_raw"].isin(control_labels),
        ],
        ["SOKM", "control"],
        default="unmapped",
    )
    output["animal"] = obs[animal_column].astype(str)
    output["known_animal"] = is_known_animal(output["animal"])
    output["state"] = obs[state_column].astype(str)
    output["sample"] = (
        obs[sample_column].astype(str)
        if sample_column in obs.columns
        else "not_available"
    )
    output["batch"] = (
        obs[batch_column].astype(str)
        if batch_column in obs.columns
        else "not_available"
    )

    if output["age_group"].eq("unmapped").any():
        bad = output.loc[output["age_group"].eq("unmapped"), "age_raw"]
        raise ValueError(f"Unmapped age labels: {bad.value_counts().to_dict()}")
    if output["analysis_role"].eq("unmapped").any():
        bad = output.loc[
            output["analysis_role"].eq("unmapped"), "treatment_raw"
        ]
        raise ValueError(f"Unmapped treatment labels: {bad.value_counts().to_dict()}")
    return output.reset_index(drop=True)
