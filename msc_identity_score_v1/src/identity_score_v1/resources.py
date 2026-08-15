"""Load and validate the frozen Identity Score v1 resources."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "gene_symbol",
    "module",
    "direction",
    "included",
}
IDENTITY_MODULE = "msc_identity_core"


def _as_bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    mapped = values.astype(str).str.strip().str.lower().map(
        {
            "true": True,
            "1": True,
            "yes": True,
            "false": False,
            "0": False,
            "no": False,
        }
    )
    if mapped.isna().any():
        raise ValueError("The `included` column must contain boolean values.")
    return mapped.astype(bool)


def read_gene_table(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    missing = REQUIRED_COLUMNS.difference(table.columns)
    if missing:
        raise ValueError(f"Gene table is missing columns: {sorted(missing)}")
    table = table.copy()
    table["gene_symbol"] = table["gene_symbol"].astype(str).str.strip()
    table["module"] = table["module"].astype(str).str.strip()
    table["direction"] = pd.to_numeric(
        table["direction"], errors="raise"
    ).astype(int)
    table["included"] = _as_bool(table["included"])
    if not set(table["direction"]).issubset({-1, 1}):
        raise ValueError("Gene directions must be -1 or +1.")
    return table


def load_frozen_resources() -> tuple[pd.DataFrame, pd.DataFrame]:
    resource_root = files("identity_score_v1").joinpath("resources")
    identity = read_gene_table(resource_root.joinpath("identity_gene_set_v1.csv"))
    all_modules = read_gene_table(
        resource_root.joinpath("reprogramming_gene_set_v1.csv")
    )

    identity_rows = identity.loc[identity["module"].eq(IDENTITY_MODULE)]
    included = identity_rows.loc[identity_rows["included"]]
    if len(included) != 11:
        raise ValueError("Identity Score v1 requires exactly 11 included genes.")
    if set(included["direction"]) != {1}:
        raise ValueError("Identity Score v1 requires 11 positive-direction genes.")
    if identity_rows.duplicated(["gene_symbol", "module"]).any():
        raise ValueError("Identity gene table contains duplicate rows.")
    return identity, all_modules
