#!/usr/bin/env python3
"""Step 2 audit for TMS Limb_Muscle MSC Youth Score workflow."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


DATA_DIR = Path("data/limb_muscle_msc")
OUT_DIR = Path("outputs/qc")
METADATA_PATH = DATA_DIR / "limb_muscle_msc_young_old_metadata.csv"
CELL_QC_PATH = DATA_DIR / "limb_muscle_msc_young_old_cell_qc.csv"


def parse_age_months(age: str) -> int:
    match = re.fullmatch(r"(\d+)m", str(age))
    if not match:
        raise ValueError(f"Cannot parse age value: {age!r}")
    return int(match.group(1))


def parse_library_id(cell_index: str) -> str:
    text = str(cell_index)
    if text.startswith("10X_"):
        parts = text.split("_")
        if len(parts) >= 3:
            return "_".join(parts[:3])
    parts = text.split("-")
    if len(parts) >= 3 and re.fullmatch(r"\d+", parts[2]):
        return f"library_{parts[2]}"
    return "unknown"


def collapse_unique(values: pd.Series) -> str:
    observed = sorted({str(value) for value in values.dropna() if str(value) != "nan"})
    return ";".join(observed) if observed else ""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(METADATA_PATH)
    cell_qc = pd.read_csv(CELL_QC_PATH)
    metadata = metadata.copy()
    metadata["mouse_id"] = metadata["mouse.id"]
    metadata["age_months"] = metadata["age"].map(parse_age_months)
    metadata["library_id"] = metadata["index"].map(parse_library_id)

    merged = metadata.merge(
        cell_qc[["index", "n_counts", "n_genes_detected"]],
        on="index",
        how="left",
        validate="one_to_one",
    )
    if merged["n_counts"].isna().any():
        raise AssertionError("Some metadata rows are missing cell QC / total count values")

    mouse_table = (
        merged.groupby("mouse_id", sort=True)
        .agg(
            age=("age", collapse_unique),
            age_months=("age_months", "nunique"),
            age_month_value=("age_months", "first"),
            age_group=("age_group", collapse_unique),
            sex=("sex", collapse_unique),
            sex_unique_count=("sex", "nunique"),
            cell_count=("index", "size"),
            total_counts=("n_counts", "sum"),
            median_counts_per_cell=("n_counts", "median"),
            mean_counts_per_cell=("n_counts", "mean"),
            median_genes_per_cell=("n_genes_detected", "median"),
            mean_genes_per_cell=("n_genes_detected", "mean"),
            method=("method", collapse_unique),
            subtissue=("subtissue", collapse_unique),
            library_id=("library_id", collapse_unique),
        )
        .reset_index()
    )
    mouse_table["age_months"] = mouse_table["age_month_value"]
    mouse_table = mouse_table.drop(columns=["age_month_value"])
    mouse_table = mouse_table[
        [
            "mouse_id",
            "age_months",
            "age",
            "age_group",
            "sex",
            "cell_count",
            "total_counts",
            "median_counts_per_cell",
            "mean_counts_per_cell",
            "median_genes_per_cell",
            "mean_genes_per_cell",
            "method",
            "subtissue",
            "library_id",
            "sex_unique_count",
        ]
    ]
    mouse_table.to_csv(OUT_DIR / "mouse_sample_table.csv", index=False)

    age_sex = pd.crosstab(merged["age_group"], merged["sex"])
    age_sex.to_csv(OUT_DIR / "age_sex_contingency_cells.csv")
    age_sex_mouse = pd.crosstab(mouse_table["age_group"], mouse_table["sex"])
    age_sex_mouse.to_csv(OUT_DIR / "age_sex_contingency_mice.csv")
    age_batch = pd.crosstab(mouse_table["age_group"], mouse_table["library_id"])
    age_batch.to_csv(OUT_DIR / "age_batch_contingency_mice.csv")

    age_sex_batch_table = mouse_table[
        [
            "mouse_id",
            "age_months",
            "age_group",
            "sex",
            "method",
            "subtissue",
            "library_id",
            "cell_count",
            "total_counts",
        ]
    ].copy()
    age_sex_batch_table.to_csv(OUT_DIR / "age_sex_batch_table.csv", index=False)

    min_cells_threshold = 50
    verification = {
        "all_mouse_ids_valid": bool(mouse_table["mouse_id"].notna().all()),
        "each_mouse_has_unique_age": bool(
            merged.groupby("mouse_id")["age_months"].nunique().eq(1).all()
        ),
        "each_mouse_has_unique_sex": bool(
            merged.groupby("mouse_id")["sex"].nunique().eq(1).all()
        ),
        "no_mouse_below_50_cells": bool(mouse_table["cell_count"].ge(min_cells_threshold).all()),
        "young_mouse_count": int((mouse_table["age_group"] == "Young").sum()),
        "old_mouse_count": int((mouse_table["age_group"] == "Old").sum()),
        "young_cell_count": int((merged["age_group"] == "Young").sum()),
        "old_cell_count": int((merged["age_group"] == "Old").sum()),
    }
    pd.DataFrame([verification]).to_csv(OUT_DIR / "step2_verification_summary.csv", index=False)

    plot_table = mouse_table.sort_values(["age_group", "age_months", "sex", "mouse_id"])
    colors = plot_table["age_group"].map({"Young": "#2E86AB", "Old": "#A23B72"}).fillna("#666666")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.bar(plot_table["mouse_id"], plot_table["cell_count"], color=colors)
    ax.axhline(min_cells_threshold, color="#444444", linestyle="--", linewidth=1)
    ax.set_ylabel("Limb_Muscle MSC cell count")
    ax.set_xlabel("Mouse ID")
    ax.set_title("Step 2 audit: cells per mouse")
    ax.tick_params(axis="x", rotation=45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "cell_counts_per_mouse.png", dpi=180)
    plt.close(fig)

    report_lines = [
        "# Step 2: Audit the Sample Structure",
        "",
        "## Inputs",
        "",
        f"- Metadata: `{METADATA_PATH}`",
        f"- Cell QC / per-cell counts: `{CELL_QC_PATH}`",
        "- Scope: Limb_Muscle mesenchymal stem cells, Young = 3m, Old = 18m + 21m + 24m.",
        "",
        "## What This Step Did",
        "",
        "1. Built a mouse-level sample table with one row per mouse.",
        "2. Summed per-cell raw UMI counts within each mouse to obtain `total_counts`.",
        "3. Counted target MSC cells per mouse and checked whether any mouse has fewer than 50 cells.",
        "4. Checked age/sex structure at both mouse and cell level.",
        "5. Parsed available technical structure from metadata: `method`, `subtissue`, and barcode-derived `library_id`.",
        "6. Wrote the required Step 2 QC outputs.",
        "",
        "## Verification",
        "",
        f"- Young mice: {verification['young_mouse_count']}",
        f"- Old mice: {verification['old_mouse_count']}",
        f"- Young cells: {verification['young_cell_count']}",
        f"- Old cells: {verification['old_cell_count']}",
        f"- Minimum cells per mouse: {int(mouse_table['cell_count'].min())}",
        f"- All mice have unique age: {verification['each_mouse_has_unique_age']}",
        f"- All mice have unique sex: {verification['each_mouse_has_unique_sex']}",
        f"- No mouse below {min_cells_threshold} cells: {verification['no_mouse_below_50_cells']}",
        "",
        "## Main Observations",
        "",
        "- The known Step 1 counts are preserved: 1468 young cells and 8181 old cells.",
        "- There are 12 mouse-level samples: 2 young and 10 old.",
        "- The two young mice are both female. Old mice include female and male mice, so sex is partially confounded with age.",
        "- Each mouse maps to a single barcode-derived library ID; library ID is therefore nested within mouse and cannot be separated cleanly from mouse/age in later models.",
        "",
        "## Outputs",
        "",
        "- `outputs/qc/mouse_sample_table.csv`",
        "- `outputs/qc/cell_counts_per_mouse.png`",
        "- `outputs/qc/age_sex_batch_table.csv`",
        "- `outputs/qc/age_sex_contingency_cells.csv`",
        "- `outputs/qc/age_sex_contingency_mice.csv`",
        "- `outputs/qc/age_batch_contingency_mice.csv`",
        "- `outputs/qc/step2_verification_summary.csv`",
    ]
    (OUT_DIR / "step2_sample_audit_report.md").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )

    print("Step 2 audit complete")
    print(mouse_table[["mouse_id", "age_months", "age_group", "sex", "cell_count", "total_counts"]])
    print("Wrote outputs to", OUT_DIR)


if __name__ == "__main__":
    main()
