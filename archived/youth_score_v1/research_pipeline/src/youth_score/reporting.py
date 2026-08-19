"""English reports and plots for trained Youth Score models."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.calibration import calibration_curve
from sklearn.metrics import roc_curve

from .config import DatasetConfig
from .data import load_prepared_bundle
from .scoring import YouthScoreEnsemble
from .utils import ensure_directory, read_json, runtime_info, utc_now, write_json


def _save_figure(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()


def generate_training_report(config: DatasetConfig) -> Path:
    output_dir = config.output_dir
    figures = ensure_directory(output_dir / "figures")
    cells = pd.read_parquet(output_dir / "oof_cell_scores.parquet")
    donors = pd.read_csv(output_dir / "oof_donor_scores.csv")
    metrics = pd.read_csv(output_dir / "model_metrics.csv")
    selection = read_json(output_dir / "selection.json")
    official = selection["official_model"]
    official_donors = donors[donors["model"].eq(official)].copy()

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(7, 5))
    sns.boxplot(data=official_donors, x="label", y="youth_score", color="#9ecae1")
    sns.stripplot(data=official_donors, x="label", y="youth_score", color="#08306b", size=7)
    plt.xticks([0, 1], ["Old", "Young"])
    plt.xlabel("Donor age group")
    plt.ylabel("Youth Score")
    plt.title(f"{config.dataset_id}: donor-level out-of-fold scores ({official})")
    _save_figure(figures / "donor_score_distribution.png")

    plt.figure(figsize=(8, 5))
    ordered = metrics.sort_values("roc_auc", ascending=False)
    sns.barplot(data=ordered, x="model", y="roc_auc", color="#4292c6")
    plt.axhline(0.5, color="black", linestyle="--", linewidth=1)
    plt.ylim(0, 1)
    plt.xlabel("Model")
    plt.ylabel("Donor-level ROC-AUC")
    plt.title(f"{config.dataset_id}: model comparison")
    plt.xticks(rotation=25, ha="right")
    _save_figure(figures / "model_auc_comparison.png")

    plt.figure(figsize=(6, 6))
    for model, group in donors.groupby("model"):
        if group["label"].nunique() != 2:
            continue
        false_positive, true_positive, _ = roc_curve(group["label"], group["youth_score"])
        auc = metrics.loc[metrics["model"].eq(model), "roc_auc"].iloc[0]
        plt.plot(false_positive, true_positive, label=f"{model} (AUC={auc:.3f})")
    plt.plot([0, 1], [0, 1], "k--", linewidth=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title(f"{config.dataset_id}: donor-level ROC")
    plt.legend(fontsize=8)
    _save_figure(figures / "donor_roc.png")

    plt.figure(figsize=(7, 5))
    sns.scatterplot(
        data=official_donors,
        x="age_months",
        y="youth_score",
        hue="label",
        size="cell_count",
        sizes=(50, 220),
        palette={0: "#ef3b2c", 1: "#2171b5"},
    )
    plt.xlabel("Chronological age (months)")
    plt.ylabel("Youth Score")
    plt.title(f"{config.dataset_id}: donor age relationship")
    _save_figure(figures / "score_vs_age.png")

    observed, predicted = calibration_curve(
        official_donors["label"],
        official_donors["youth_score"],
        n_bins=min(5, max(2, len(official_donors) // 3)),
        strategy="quantile",
    )
    plt.figure(figsize=(6, 6))
    plt.plot([0, 1], [0, 1], "k--", linewidth=1, label="Ideal")
    plt.plot(predicted, observed, marker="o", color="#2171b5", label=official)
    plt.xlabel("Mean predicted Youth Score")
    plt.ylabel("Observed young fraction")
    plt.title(f"{config.dataset_id}: donor-level calibration")
    plt.legend()
    _save_figure(figures / "donor_calibration.png")

    bundle = load_prepared_bundle(config.processed_dir)
    expected = bundle.cells.groupby("age_group").agg(cells=("cell_id", "size"), mice=("mouse_id", "nunique"))
    metric_markdown = metrics.round(4).to_markdown(index=False)
    report = f"""# Youth Score Training Report: {config.dataset_id}

Generated: {utc_now()}

## Dataset

- Role: `{config.role}`
- Modality: `{config.modality}`
- Tissue: `{config.tissue}`
- Cell type: `{config.cell_type}`
- Young ages: `{', '.join(config.young_ages)}`
- Old ages: `{', '.join(config.old_ages)}`
- Cells after matrix-derived QC: `{len(bundle.cells)}`
- Donors: `{bundle.cells['mouse_id'].nunique()}`
- Young cells/donors: `{int(expected.loc['young', 'cells'])}` / `{int(expected.loc['young', 'mice'])}`
- Old cells/donors: `{int(expected.loc['old', 'cells'])}` / `{int(expected.loc['old', 'mice'])}`

## Model Selection

- Official model: `{official}`
- Evidence status: `{selection['status']}`
- Rule: {selection['selection_rule']}
- Technical-only AUC: `{selection['sensitivity']['technical_auc']:.4f}`
- Official-model AUC: `{selection['sensitivity']['official_auc']:.4f}`

## Donor-level Metrics

{metric_markdown}

## Interpretation

The score is a cross-validated estimate of similarity to the TMS 3-month reference state. A higher score means more young-like expression within this tissue/cell-type model. It is not a safety score, causal rejuvenation measurement, clinical age, or treatment recommendation.

The primary unit of validation is the mouse, not the cell. Sex, sequencing plate, library complexity, and age are partially confounded in TMS. The `confound_limited` status means the model should remain exploratory even if its discrimination metric is high.

## Files

- `oof_cell_scores.parquet`: out-of-fold cell scores.
- `oof_donor_scores.csv`: donor-aggregated scores.
- `model_metrics.csv`: baseline and Transformer comparison.
- `selection.json`: selected model and confounding status.
- `folds/`: fold-specific tokenizers, calibration, baselines, and Transformer checkpoints.
- `figures/`: report plots.
"""
    report_path = output_dir / "training_report.md"
    report_path.write_text(report, encoding="utf-8")
    write_json(
        output_dir / "training_report.json",
        {
            "created_at": utc_now(),
            "dataset_id": config.dataset_id,
            "cohort": {
                "cells": len(bundle.cells),
                "donors": int(bundle.cells["mouse_id"].nunique()),
                "young_cells": int(expected.loc["young", "cells"]),
                "young_donors": int(expected.loc["young", "mice"]),
                "old_cells": int(expected.loc["old", "cells"]),
                "old_donors": int(expected.loc["old", "mice"]),
            },
            "selection": selection,
            "metrics": metrics.to_dict(orient="records"),
            "figures": sorted(path.name for path in figures.glob("*.png")),
        },
    )
    return report_path


def score_limb_droplet_sensitivity(
    sensitivity_config: DatasetConfig,
    trained_model_directory: str | Path,
    device: str = "auto",
) -> Path:
    bundle = load_prepared_bundle(sensitivity_config.processed_dir)
    scorer = YouthScoreEnsemble(trained_model_directory, device=device)
    scores = scorer.score_counts(
        bundle.counts,
        bundle.genes["gene_name"].astype(str).tolist(),
        bundle.cells["cell_id"].astype(str).tolist(),
    )
    scores = scores.merge(
        bundle.cells[["cell_id", "mouse_id", "age", "age_months", "sex"]],
        on="cell_id",
        how="left",
        validate="one_to_one",
    )
    output_dir = ensure_directory(
        sensitivity_config.project_root / "outputs" / "youth_score" / sensitivity_config.dataset_id
    )
    scores.to_parquet(output_dir / "cell_scores.parquet", index=False)
    donor = (
        scores.groupby(["mouse_id", "age", "age_months", "sex"], as_index=False)
        .agg(youth_score=("youth_score", "mean"), predicted_age_months=("predicted_age_months", "mean"), cells=("cell_id", "size"))
        .sort_values(["age_months", "mouse_id"])
    )
    donor.to_csv(output_dir / "donor_scores.csv", index=False)
    age_summary = donor.groupby(["age", "age_months"], as_index=False).agg(
        youth_score=("youth_score", "mean"), donors=("mouse_id", "nunique"), cells=("cells", "sum")
    )
    age_summary.to_csv(output_dir / "age_summary.csv", index=False)
    age_correlation = float(spearmanr(donor["age_months"], donor["youth_score"]).statistic)

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 5))
    sns.lineplot(data=donor, x="age_months", y="youth_score", marker="o", errorbar=None)
    sns.scatterplot(data=donor, x="age_months", y="youth_score", hue="sex", size="cells", sizes=(50, 220))
    plt.xlabel("Chronological age (months)")
    plt.ylabel("Youth Score")
    plt.title("Droplet Limb Muscle MSC cross-modality sensitivity")
    _save_figure(output_dir / "score_by_age.png")
    write_json(
        output_dir / "sensitivity_summary.json",
        {
            "created_at": utc_now(),
            "source_model": str(trained_model_directory),
            "official_model": scorer.official_model,
            "minimum_gene_overlap": float(scores["gene_overlap"].min()),
            "donor_age_spearman": age_correlation,
            "age_direction_pass": bool(np.isfinite(age_correlation) and age_correlation < 0),
            "age_summary": age_summary.to_dict(orient="records"),
            "interpretation": "Sensitivity analysis only; age, sex, donor count, and modality are confounded.",
        },
    )
    return output_dir


def generate_project_summary(project_root: str | Path) -> Path:
    project_root = Path(project_root).resolve()
    output_root = ensure_directory(project_root / "outputs" / "youth_score")
    internal = {}
    for dataset_id in ("scat_facs", "limb_facs"):
        directory = output_root / dataset_id
        if not (directory / "selection.json").exists():
            continue
        internal[dataset_id] = {
            "selection": read_json(directory / "selection.json"),
            "metrics": pd.read_csv(directory / "model_metrics.csv").to_dict(orient="records"),
        }
        all_genes = directory / "all_genes_sensitivity"
        if (all_genes / "selection.json").exists():
            internal[dataset_id]["all_genes_sensitivity"] = {
                "selection": read_json(all_genes / "selection.json"),
                "metrics": pd.read_csv(all_genes / "model_metrics.csv").to_dict(orient="records"),
            }
    sensitivity_path = output_root / "limb_droplet_sensitivity" / "sensitivity_summary.json"
    external_path = output_root / "external" / "GSE176206" / "combined_summary.json"
    report = {
        "created_at": utc_now(),
        "runtime": runtime_info(),
        "internal_validation": internal,
        "droplet_sensitivity": read_json(sensitivity_path) if sensitivity_path.exists() else None,
        "external_validation": read_json(external_path) if external_path.exists() else None,
        "scientific_scope": (
            "Research-only transcriptomic similarity score; not a safety score, causal "
            "rejuvenation measurement, clinical age, or treatment recommendation."
        ),
    }
    target = output_root / "run_report.json"
    write_json(target, report)
    return target
