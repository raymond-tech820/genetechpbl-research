from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def youth_validation_figure() -> None:
    metrics = pd.read_csv(DATA / "model_metrics.csv")
    droplet = pd.read_csv(DATA / "tms_droplet_donor_scores.csv")

    label_map = {
        "gene_signature": "Gene\nsignature",
        "elastic_net": "Elastic net",
        "gene_transformer": "Gene-token\nTransformer",
        "technical_only": "Technical-only\ndiagnostic",
    }
    colors = ["#1F77B4", "#4C78A8", "#72B7B2", "#B8B8B8"]

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.65), gridspec_kw={"wspace": 0.30})

    ax = axes[0]
    x = np.arange(len(metrics))
    y = metrics["roc_auc"].to_numpy()
    low = y - metrics["roc_auc_ci_low"].to_numpy()
    high = metrics["roc_auc_ci_high"].to_numpy() - y
    bars = ax.bar(x, y, color=colors, width=0.68, zorder=2)
    ax.errorbar(x, y, yerr=np.vstack([low, high]), fmt="none", ecolor="#303030", capsize=3, lw=1, zorder=3)
    ax.axhline(0.5, color="#777777", ls="--", lw=1)
    ax.set_ylim(0.45, 1.035)
    ax.set_ylabel("Donor-level out-of-fold ROC-AUC")
    ax.set_xticks(x, [label_map[m] for m in metrics["model"]])
    ax.set_title("A  Donor-grouped model evaluation", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#E6E6E6", lw=0.8, zorder=0)
    for bar, value in zip(bars, y):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.017, f"{value:.3f}", ha="center", va="bottom", fontsize=8)

    ax = axes[1]
    age_palette = {
        1.0: "#2CA02C",
        3.0: "#66A61E",
        18.0: "#E6AB02",
        21.0: "#E67E22",
        24.0: "#D95F02",
        30.0: "#A63603",
    }
    for age, frame in droplet.groupby("age_months", sort=True):
        jitter = np.linspace(-0.22, 0.22, len(frame)) if len(frame) > 1 else np.array([0.0])
        ax.scatter(
            frame["age_months"].to_numpy() + jitter,
            frame["youth_score"],
            s=43,
            color=age_palette[float(age)],
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
        ax.plot(
            [age - 0.36, age + 0.36],
            [frame["youth_score"].mean()] * 2,
            color="#222222",
            lw=1.6,
            zorder=4,
        )
    rho, _ = spearmanr(droplet["age_months"], droplet["youth_score"])
    ax.set_xlabel("Age (months)")
    ax.set_ylabel("Youth Score")
    ax.set_xticks(sorted(droplet["age_months"].unique()))
    ax.set_ylim(0.0, 0.60)
    ax.set_title("B  Within-TMS droplet sensitivity", loc="left", fontweight="bold")
    ax.text(0.97, 0.95, f"Spearman rho = {rho:.3f}\nn = {droplet['mouse_id'].nunique()} animals", transform=ax.transAxes, ha="right", va="top")
    ax.grid(axis="y", color="#E6E6E6", lw=0.8, zorder=0)

    fig.savefig(FIGURES / "youth_validation_and_cross_assay.pdf")
    fig.savefig(FIGURES / "youth_validation_and_cross_assay.png", dpi=300)
    plt.close(fig)


def identity_exact_arm_figure() -> None:
    data = pd.read_csv(DATA / "identity_donor_by_condition.csv")
    data["arm_label"] = data["exact_treatment_arm"].map(
        {
            "Tg-/Dox+": "Tg-/Dox+\ncontrol",
            "Tg+/Dox-": "Tg+/Dox-\ncontrol",
            "Tg+/Dox+": "Tg+/Dox+\nSOKM",
        }
    )
    arm_order = ["Tg-/Dox+", "Tg+/Dox-", "Tg+/Dox+"]
    age_order = ["young", "aged"]
    arm_colors = {"Tg-/Dox+": "#7A9CC6", "Tg+/Dox-": "#63B7AF", "Tg+/Dox+": "#D95F5F"}
    offsets = {"young": -0.19, "aged": 0.19}
    age_markers = {"young": "o", "aged": "s"}

    panels = [
        ("identity_primary_median", "Background-corrected primary score", "A  Primary Identity Score"),
        ("identity_rank_median", "Rank-based sensitivity score", "B  Rank-based Identity Score"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.85), gridspec_kw={"wspace": 0.28})

    for ax, (value_col, y_label, title) in zip(axes, panels):
        for arm_index, arm in enumerate(arm_order):
            for age in age_order:
                frame = data[(data["exact_treatment_arm"] == arm) & (data["age_group"] == age)].copy()
                values = frame[value_col].to_numpy()
                x_center = arm_index + offsets[age]
                jitter = np.linspace(-0.045, 0.045, len(values))
                ax.scatter(
                    x_center + jitter,
                    values,
                    color=arm_colors[arm],
                    marker=age_markers[age],
                    s=47,
                    edgecolor="white",
                    linewidth=0.75,
                    zorder=3,
                    label=age.capitalize() if arm_index == 0 else None,
                )
                median = float(np.median(values))
                ax.plot([x_center - 0.105, x_center + 0.105], [median, median], color="#1F1F1F", lw=1.8, zorder=4)
        ax.set_xticks(range(len(arm_order)), [data.loc[data["exact_treatment_arm"] == arm, "arm_label"].iloc[0] for arm in arm_order])
        ax.set_ylabel(y_label)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.grid(axis="y", color="#E6E6E6", lw=0.8, zorder=0)
        ax.margins(x=0.10)

    axes[0].set_ylim(0.08, 0.33)
    axes[1].set_ylim(0.27, 0.54)
    handles, labels = axes[1].get_legend_handles_labels()
    axes[1].legend(handles, labels, title="Age group", frameon=False, loc="lower left")
    fig.text(0.5, -0.025, "Points represent animal-level medians (n = 3 animals per exact arm); horizontal bars show the condition median.", ha="center", fontsize=8.6)
    fig.savefig(FIGURES / "identity_exact_arm_scores.pdf")
    fig.savefig(FIGURES / "identity_exact_arm_scores.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    set_style()
    youth_validation_figure()
    identity_exact_arm_figure()
