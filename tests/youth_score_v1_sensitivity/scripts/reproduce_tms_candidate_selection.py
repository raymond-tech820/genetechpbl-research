#!/usr/bin/env python3
"""Reproduce TMS tissue-cell candidate rankings from metadata.

The candidate selection report ranks tissue/cell-type groups with:

Young = 3m
Old   = 18m + 21m + 24m

The expression matrix is not needed for the ranking itself, but this script can
check that a BPCells directory has the expected column names for the metadata.
"""

from __future__ import annotations

import argparse
import math
import struct
from pathlib import Path

import numpy as np
import pandas as pd


YOUNG_AGES = {"3m"}
OLD_AGES = {"18m", "21m", "24m"}
GROUP_COLUMNS = ["tissue", "cell_ontology_class"]
DISPLAY_COLUMNS = [
    "rank",
    "tissue",
    "cell_ontology_class",
    "young_cells",
    "old_cells",
    "min_group_size",
    "young_old_balance",
    "young_mice",
    "old_mice",
    "score",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/tabula-muris-senis-droplet_obs.csv"),
        help="TMS metadata CSV with age, tissue, cell_ontology_class, and mouse.id columns.",
    )
    parser.add_argument(
        "--bpcells-dir",
        type=Path,
        default=Path("data/tabula-muris-senis-droplet_bpcells"),
        help="Optional BPCells matrix directory used for shape/name consistency checks.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("outputs/tms_candidate_selection"),
        help="Directory for generated ranking CSV and Markdown files.",
    )
    parser.add_argument(
        "--prefix",
        default="droplet",
        help="Filename prefix for outputs, e.g. droplet or facs.",
    )
    return parser.parse_args()


def read_bpcells_shape(shape_path: Path) -> tuple[int, int]:
    data = shape_path.read_bytes()
    magic = b"UINT32v1"
    if not data.startswith(magic):
        raise ValueError(f"Unsupported BPCells shape header in {shape_path}")
    rows, cols = struct.unpack("<II", data[len(magic) : len(magic) + 8])
    return rows, cols


def validate_bpcells(metadata: pd.DataFrame, bpcells_dir: Path) -> dict[str, object]:
    if not bpcells_dir.exists():
        return {"available": False, "message": f"{bpcells_dir} does not exist"}

    shape = read_bpcells_shape(bpcells_dir / "shape")
    row_count = sum(1 for _ in (bpcells_dir / "row_names").open())
    col_names = pd.read_csv(bpcells_dir / "col_names", header=None, names=["index"])

    metadata_index = metadata["index"].astype(str).reset_index(drop=True)
    col_index = col_names["index"].astype(str)
    same_order = bool(metadata_index.equals(col_index))

    return {
        "available": True,
        "shape_rows": shape[0],
        "shape_cols": shape[1],
        "row_names": row_count,
        "col_names": len(col_names),
        "metadata_rows": len(metadata),
        "metadata_matches_col_names": same_order,
    }


def build_candidate_table(metadata: pd.DataFrame) -> pd.DataFrame:
    required = {"age", "tissue", "cell_ontology_class", "mouse.id"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing required columns: {sorted(missing)}")

    scoped = metadata.loc[metadata["age"].isin(YOUNG_AGES | OLD_AGES)].copy()
    scoped["age_group"] = np.where(scoped["age"].isin(YOUNG_AGES), "Young", "Old")

    cell_counts = (
        scoped.groupby(GROUP_COLUMNS + ["age_group"], dropna=False)
        .size()
        .unstack("age_group", fill_value=0)
        .rename(columns={"Young": "young_cells", "Old": "old_cells"})
    )

    for column in ["young_cells", "old_cells"]:
        if column not in cell_counts:
            cell_counts[column] = 0

    mouse_counts = (
        scoped.groupby(GROUP_COLUMNS + ["age_group"], dropna=False)["mouse.id"]
        .nunique()
        .unstack("age_group", fill_value=0)
        .rename(columns={"Young": "young_mice", "Old": "old_mice"})
    )

    for column in ["young_mice", "old_mice"]:
        if column not in mouse_counts:
            mouse_counts[column] = 0

    table = cell_counts.join(mouse_counts).reset_index()
    table = table[(table["young_cells"] > 0) & (table["old_cells"] > 0)].copy()

    table["min_group_size"] = table[["young_cells", "old_cells"]].min(axis=1)
    table["young_old_balance"] = (
        2 * table["min_group_size"] / (table["young_cells"] + table["old_cells"])
    )
    table["min_mice"] = table[["young_mice", "old_mice"]].min(axis=1)
    table["mouse_balance"] = (
        2 * table["min_mice"] / (table["young_mice"] + table["old_mice"])
    )

    max_log_min = math.log1p(table["min_group_size"].max())
    max_log_mice = math.log1p(table["min_mice"].max())
    table["normalized_sample_size"] = np.log1p(table["min_group_size"]) / max_log_min
    table["normalized_mouse_replicates"] = np.log1p(table["min_mice"]) / max_log_mice
    table["composite_score"] = 100 * (
        0.50 * table["normalized_sample_size"]
        + 0.25 * table["young_old_balance"]
        + 0.25 * table["normalized_mouse_replicates"]
    )

    table = table.sort_values(
        ["min_group_size", "young_old_balance", "tissue", "cell_ontology_class"],
        ascending=[False, False, True, True],
    )
    return table.reset_index(drop=True)


def ranked_view(table: pd.DataFrame, score_column: str) -> pd.DataFrame:
    view = table.copy()
    view["score"] = view[score_column]
    view = view.sort_values(
        [
            "score",
            "min_group_size",
            "old_cells",
            "young_old_balance",
            "tissue",
            "cell_ontology_class",
        ],
        ascending=[False, False, False, False, True, True],
    ).reset_index(drop=True)
    view.insert(0, "rank", np.arange(1, len(view) + 1))
    return view[DISPLAY_COLUMNS]


def dataframe_to_markdown(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    rows = [[str(value) for value in row] for row in frame.itertuples(index=False, name=None)]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def write_markdown_report(
    out_path: Path, prefix: str, rankings: dict[str, pd.DataFrame], validation: dict[str, object]
) -> None:
    sections = [
        f"# TMS candidate selection reproduction: {prefix}",
        "",
        "Age rule: Young = 3m; Old = 18m + 21m + 24m.",
        "",
        "## BPCells consistency check",
        "",
    ]

    if validation.get("available"):
        sections.extend(
            [
                f"- matrix shape: {validation['shape_rows']} genes x {validation['shape_cols']} cells",
                f"- row_names: {validation['row_names']}",
                f"- col_names: {validation['col_names']}",
                f"- metadata rows: {validation['metadata_rows']}",
                f"- metadata index matches col_names order: {validation['metadata_matches_col_names']}",
                "",
            ]
        )
    else:
        sections.extend([f"- {validation['message']}", ""])

    labels = {
        "minimum_group_size": "Top by minimum group size",
        "cell_balance": "Top by cell balance",
        "donor_replication": "Top by donor replication",
        "composite": "Top composite score",
    }
    rename = {
        "rank": "Rank",
        "tissue": "Tissue",
        "cell_ontology_class": "Cell type",
        "young_cells": "Young",
        "old_cells": "Old",
        "min_group_size": "Min",
        "young_old_balance": "Young/Old balance",
        "young_mice": "Young mice",
        "old_mice": "Old mice",
        "score": "Score",
    }
    for key, label in labels.items():
        top = rankings[key].head(10).copy()
        top["young_old_balance"] = top["young_old_balance"].round(3)
        top["score"] = top["score"].round(2)
        sections.extend(
            ["## " + label, "", dataframe_to_markdown(top.rename(columns=rename)), ""]
        )

    out_path.write_text("\n".join(sections), encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata, low_memory=False)
    validation = validate_bpcells(metadata, args.bpcells_dir)
    candidates = build_candidate_table(metadata)

    candidates.to_csv(args.out_dir / f"{args.prefix}_candidate_scores.csv", index=False)

    scored = candidates.assign(
        cell_balance_score=candidates["min_group_size"] * candidates["young_old_balance"],
        donor_replication_score=candidates["min_group_size"] * candidates["mouse_balance"],
    )
    rankings = {
        "minimum_group_size": ranked_view(scored, "min_group_size"),
        "cell_balance": ranked_view(scored, "cell_balance_score"),
        "donor_replication": ranked_view(scored, "donor_replication_score"),
        "composite": ranked_view(scored, "composite_score"),
    }

    for key, ranking in rankings.items():
        ranking.to_csv(args.out_dir / f"{args.prefix}_{key}_top.csv", index=False)

    write_markdown_report(
        args.out_dir / f"{args.prefix}_candidate_selection_reproduction.md",
        args.prefix,
        rankings,
        validation,
    )

    print(f"Read {len(metadata):,} metadata rows")
    if validation.get("available"):
        print(
            "BPCells:",
            f"{validation['shape_rows']:,} genes x {validation['shape_cols']:,} cells;",
            "metadata order match =",
            validation["metadata_matches_col_names"],
        )
    print(f"Wrote outputs to {args.out_dir}")


if __name__ == "__main__":
    main()
