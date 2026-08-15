"""Generate publication-style figures for the Youth Score progress summary."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = ROOT / "outputs" / "youth_score"
OUTPUT_DIR = RESULT_ROOT / "progress_summary" / "figures"

MODEL_ORDER = ["gene_signature", "elastic_net", "gene_transformer", "technical_only"]
MODEL_LABELS = {
    "gene_signature": "Gene signature",
    "elastic_net": "Elastic Net",
    "gene_transformer": "Transformer",
    "technical_only": "Technical only",
}
MODEL_COLORS = {
    "gene_signature": "#009E73",
    "elastic_net": "#0072B2",
    "gene_transformer": "#CC79A7",
    "technical_only": "#7A7A7A",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT_DIR / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_internal_auc() -> None:
    datasets = [
        ("SCAT adipose MSC", "scat_facs", "elastic_net", 15),
        ("Limb Muscle MSC", "limb_facs", "gene_signature", 14),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.8), sharex=True, sharey=True)
    y_positions = np.arange(len(MODEL_ORDER))[::-1]

    for axis, (title, dataset_id, selected, donor_count) in zip(axes, datasets):
        metrics = pd.read_csv(RESULT_ROOT / dataset_id / "model_metrics.csv").set_index("model")
        for y, model in zip(y_positions, MODEL_ORDER):
            row = metrics.loc[model]
            value = float(row["roc_auc"])
            lower = float(row["roc_auc_ci_low"])
            upper = float(row["roc_auc_ci_high"])
            axis.errorbar(
                value,
                y,
                xerr=np.array([[value - lower], [upper - value]]),
                fmt="o",
                markersize=8 if model == selected else 6,
                markerfacecolor=MODEL_COLORS[model],
                markeredgecolor="black" if model == selected else "white",
                markeredgewidth=1.2 if model == selected else 0.7,
                ecolor=MODEL_COLORS[model],
                elinewidth=1.6,
                capsize=3,
                zorder=3,
            )
            axis.text(min(1.015, value + 0.018), y + 0.14, f"{value:.3f}", fontsize=8, ha="left")

        axis.axvline(0.5, color="#555555", linestyle="--", linewidth=1, zorder=1)
        axis.set_title(f"{title}\n(n = {donor_count} donors)")
        axis.set_xlim(0.30, 1.045)
        axis.set_xticks(np.arange(0.4, 1.01, 0.1))
        axis.grid(axis="x", color="#E6E6E6", linewidth=0.8, zorder=0)
        axis.set_xlabel("Donor-level ROC-AUC")

    axes[0].set_yticks(y_positions, [MODEL_LABELS[model] for model in MODEL_ORDER])
    for tick, model in zip(axes[0].get_yticklabels(), MODEL_ORDER):
        if model == "elastic_net":
            tick.set_fontweight("bold")
    axes[1].tick_params(axis="y", labelleft=False)
    fig.suptitle("Internal validation by donor-grouped five-fold cross-validation", y=1.03, fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        -0.01,
        "Points are pooled out-of-fold donor ROC-AUC values; bars are 95% donor-bootstrap intervals. "
        "Black-edged points denote the selected model; dashed line denotes chance.",
        ha="center",
        va="top",
        fontsize=8.5,
    )
    fig.tight_layout()
    save_figure(fig, "figure_1_internal_auc")


def plot_droplet_age_trajectory() -> None:
    donor = pd.read_csv(RESULT_ROOT / "limb_droplet_sensitivity" / "donor_scores.csv")
    with (RESULT_ROOT / "limb_droplet_sensitivity" / "sensitivity_summary.json").open(encoding="utf-8") as handle:
        summary = json.load(handle)

    ages = sorted(donor["age_months"].unique())
    age_means = donor.groupby("age_months", sort=True)["youth_score"].mean()
    fig, axis = plt.subplots(figsize=(7.8, 4.8))

    for age in ages:
        values = donor.loc[donor["age_months"].eq(age), "youth_score"].to_numpy()
        offsets = np.linspace(-0.35, 0.35, len(values)) if len(values) > 1 else np.array([0.0])
        axis.scatter(
            age + offsets,
            values,
            s=44,
            facecolor="#56B4E9",
            edgecolor="white",
            linewidth=0.8,
            alpha=0.95,
            zorder=3,
        )
        axis.text(age, 0.018, f"n={len(values)}", ha="center", va="bottom", fontsize=8, color="#444444")

    axis.plot(
        age_means.index,
        age_means.values,
        color="#0072B2",
        linewidth=2.0,
        marker="o",
        markersize=6,
        markerfacecolor="white",
        markeredgewidth=1.5,
        zorder=4,
        label="Age-group mean",
    )
    axis.set_xticks(ages, [f"{int(age)}m" for age in ages])
    axis.set_ylim(0.0, 0.56)
    axis.set_xlabel("Mouse age")
    axis.set_ylabel("Youth Score")
    axis.set_title("Cross-modality sensitivity: FACS Limb model applied to Droplet Limb MSC")
    axis.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    axis.text(
        0.98,
        0.95,
        f"Donor-level Spearman ρ = {summary['donor_age_spearman']:.3f}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#BDBDBD"},
    )
    axis.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.84))
    fig.text(
        0.5,
        -0.01,
        "Small points are individual donors; the line connects unweighted donor means. "
        "This is a sensitivity analysis: age, sex, donor count, and assay modality are partially confounded.",
        ha="center",
        va="top",
        fontsize=8.5,
    )
    fig.tight_layout()
    save_figure(fig, "figure_2_droplet_age_trajectory")


def load_external() -> tuple[dict[str, pd.DataFrame], dict]:
    replicates = {
        "Adipo": pd.read_csv(RESULT_ROOT / "external" / "GSE176206" / "adipo_sokm" / "replicate_scores.csv"),
        "MSC": pd.read_csv(RESULT_ROOT / "external" / "GSE176206" / "msc_sokm" / "replicate_scores.csv"),
    }
    with (RESULT_ROOT / "external" / "GSE176206" / "combined_summary.json").open(encoding="utf-8") as handle:
        summary = json.load(handle)
    return replicates, summary


def scatter_condition(
    axis: plt.Axes,
    dataset_index: int,
    values: np.ndarray,
    offset: float,
    color: str,
    marker: str,
) -> None:
    jitter = np.linspace(-0.035, 0.035, len(values)) if len(values) > 1 else np.array([0.0])
    x_values = dataset_index + offset + jitter
    axis.scatter(
        x_values,
        values,
        s=48,
        marker=marker,
        facecolor=color if len(values) > 1 else "white",
        edgecolor=color,
        linewidth=1.2,
        zorder=3,
    )
    mean_value = float(np.mean(values))
    axis.plot(
        [dataset_index + offset - 0.09, dataset_index + offset + 0.09],
        [mean_value, mean_value],
        color=color,
        linewidth=2.2,
        zorder=4,
    )


def plot_external_age_validation() -> None:
    replicates, summary = load_external()
    datasets = ["Adipo", "MSC"]
    colors = {"young": "#0072B2", "aged": "#D55E00"}
    offsets = {"young": -0.16, "aged": 0.16}
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), gridspec_kw={"width_ratios": [1.35, 1.0]})

    for dataset_index, dataset in enumerate(datasets):
        frame = replicates[dataset]
        for age_group in ("young", "aged"):
            values = frame.loc[
                frame["external_age_group"].eq(age_group) & frame["external_treatment"].eq("control"),
                "youth_score",
            ].to_numpy()
            scatter_condition(
                axes[0],
                dataset_index,
                values,
                offsets[age_group],
                colors[age_group],
                "o" if len(values) > 1 else "D",
            )

    axes[0].set_xticks(range(len(datasets)), datasets)
    axes[0].set_xlim(-0.5, 1.5)
    axes[0].set_ylim(0.40, 0.52)
    axes[0].set_ylabel("Youth Score")
    axes[0].set_title("Control-condition scores")
    axes[0].grid(axis="y", color="#E6E6E6", linewidth=0.8)
    axes[0].legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["young"], markeredgecolor=colors["young"], label="Young control"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["aged"], markeredgecolor=colors["aged"], label="Aged control"),
        ],
        frameon=False,
        loc="lower left",
    )

    contrast_rows = [
        (
            "Adipo",
            summary["adipo_sokm"]["contrasts"]["young_control_minus_aged_control"],
            summary["adipo_sokm"]["contrast_intervals_95"]["young_control_minus_aged_control"],
        ),
        (
            "MSC",
            summary["msc_sokm"]["contrasts"]["young_control_minus_aged_control"],
            summary["msc_sokm"]["contrast_intervals_95"]["young_control_minus_aged_control"],
        ),
    ]
    for y, (dataset, value, interval) in zip([1, 0], contrast_rows):
        if interval is None:
            axes[1].scatter(value, y, marker="D", s=58, facecolor="white", edgecolor="#222222", linewidth=1.3, zorder=3)
            axes[1].text(value + 0.002, y + 0.12, "No donor CI", fontsize=8, ha="left")
        else:
            axes[1].errorbar(
                value,
                y,
                xerr=np.array([[value - interval[0]], [interval[1] - value]]),
                fmt="o",
                color="#0072B2",
                markersize=7,
                capsize=4,
                linewidth=1.7,
                zorder=3,
            )
        label_y = y - 0.18 if y == 1 else y + 0.13
        axes[1].text(value, label_y, f"{value:+.3f}", ha="center", fontsize=8.5)

    axes[1].axvline(0.0, color="#555555", linestyle="--", linewidth=1)
    axes[1].set_yticks([1, 0], ["Adipo", "MSC"])
    axes[1].set_xlim(-0.008, 0.063)
    axes[1].set_xlabel("Young control − Aged control")
    axes[1].set_title("External age contrast")
    axes[1].grid(axis="x", color="#E6E6E6", linewidth=0.8)
    fig.suptitle("External age-direction validation in GSE176206", y=1.03, fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        -0.01,
        "MSC points represent three animals per condition and the contrast bar is a 95% donor-bootstrap interval. "
        "Adipo has one pooled library per condition in the processed metadata, so no donor-level interval is available.",
        ha="center",
        va="top",
        fontsize=8.5,
    )
    fig.tight_layout()
    save_figure(fig, "figure_3_external_age_validation")


def plot_external_sokm_contrast() -> None:
    replicates, summary = load_external()
    datasets = ["Adipo", "MSC"]
    colors = {"control": "#4D4D4D", "SOKM": "#CC79A7"}
    offsets = {"control": -0.16, "SOKM": 0.16}
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), gridspec_kw={"width_ratios": [1.35, 1.0]})

    for dataset_index, dataset in enumerate(datasets):
        frame = replicates[dataset]
        for treatment in ("control", "SOKM"):
            values = frame.loc[
                frame["external_age_group"].eq("aged") & frame["external_treatment"].eq(treatment),
                "youth_score",
            ].to_numpy()
            scatter_condition(
                axes[0],
                dataset_index,
                values,
                offsets[treatment],
                colors[treatment],
                "o" if len(values) > 1 else "D",
            )

    axes[0].set_xticks(range(len(datasets)), datasets)
    axes[0].set_xlim(-0.5, 1.5)
    axes[0].set_ylim(0.27, 0.49)
    axes[0].set_ylabel("Youth Score")
    axes[0].set_title("Aged-condition scores")
    axes[0].grid(axis="y", color="#E6E6E6", linewidth=0.8)
    axes[0].legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["control"], markeredgecolor=colors["control"], label="Aged control"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["SOKM"], markeredgecolor=colors["SOKM"], label="Aged SOKM"),
        ],
        frameon=False,
        loc="lower left",
    )

    contrast_rows = [
        (
            "Adipo",
            summary["adipo_sokm"]["contrasts"]["aged_sokm_minus_aged_control"],
            summary["adipo_sokm"]["contrast_intervals_95"]["aged_sokm_minus_aged_control"],
        ),
        (
            "MSC",
            summary["msc_sokm"]["contrasts"]["aged_sokm_minus_aged_control"],
            summary["msc_sokm"]["contrast_intervals_95"]["aged_sokm_minus_aged_control"],
        ),
    ]
    for y, (dataset, value, interval) in zip([1, 0], contrast_rows):
        if interval is None:
            axes[1].scatter(value, y, marker="D", s=58, facecolor="white", edgecolor="#222222", linewidth=1.3, zorder=3)
            axes[1].text(value + 0.006, y + 0.12, "No donor CI", fontsize=8, ha="left")
        else:
            axes[1].errorbar(
                value,
                y,
                xerr=np.array([[value - interval[0]], [interval[1] - value]]),
                fmt="o",
                color="#CC79A7",
                markersize=7,
                capsize=4,
                linewidth=1.7,
                zorder=3,
            )
        label_y = y - 0.18 if y == 1 else y + 0.13
        axes[1].text(value, label_y, f"{value:+.3f}", ha="center", fontsize=8.5)

    axes[1].axvline(0.0, color="#555555", linestyle="--", linewidth=1)
    axes[1].set_yticks([1, 0], ["Adipo", "MSC"])
    axes[1].set_xlim(-0.185, 0.035)
    axes[1].set_xticks([-0.15, -0.10, -0.05, 0.00])
    axes[1].set_xlabel("Aged SOKM − Aged control")
    axes[1].set_title("SOKM effect contrast")
    axes[1].grid(axis="x", color="#E6E6E6", linewidth=0.8)
    fig.suptitle("External SOKM analysis in GSE176206", y=1.03, fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        -0.01,
        "Positive values would support movement toward the TMS-defined youth axis. "
        "The MSC bar is a 95% donor-bootstrap interval; Adipo lacks donor-resolved replication.",
        ha="center",
        va="top",
        fontsize=8.5,
    )
    fig.tight_layout()
    save_figure(fig, "figure_4_external_sokm_contrast")


def main() -> None:
    configure_style()
    plot_internal_auc()
    plot_droplet_age_trajectory()
    plot_external_age_validation()
    plot_external_sokm_contrast()
    print(f"Figures written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
