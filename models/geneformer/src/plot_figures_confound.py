"""
Publication figures from existing axis-validation results. Plotting only, no new analysis.

Figure 1: per-cell scatter of axis_position vs n_genes_detected, coloured by age group,
           from results/axis_validation_per_cell_axis_position.csv.
Figure 2: forest-style summary of donor LOO accuracy and permutation p across the three
           section 6.6 conditions (raw / n_genes-regressed-out / matched depth band), from
           results/axis_validation_summary.json and results/axis_validation_deconfound.json.

Run:    python src/plot_figures_confound.py
Output: results/figures/fig1_axis_vs_ngenes_confound.{pdf,png}
        results/figures/fig2_deconfound_forest.{pdf,png}
"""

import json

import matplotlib.pyplot as plt
import pandas as pd

RESULTS_DIR = "results"
FIG_DIR = "results/figures"
DPI = 300

# Colourblind-validated categorical pair (skill palette slots 1-2, adjacent-pair CVD-safe)
COLOR_OLD = "#2a78d6"    # blue
COLOR_YOUNG = "#eb6834"  # orange
COLOR_INK = "#0b0b0b"
COLOR_MUTED = "#52514e"
COLOR_GRID = "#d8d7d2"
COLOR_REF = "#e34948"    # red, reserved for the p=0.05 reference line only

plt.rcParams.update({
    "font.size": 11,
    "axes.edgecolor": COLOR_MUTED,
    "axes.labelcolor": COLOR_INK,
    "text.color": COLOR_INK,
    "xtick.color": COLOR_MUTED,
    "ytick.color": COLOR_MUTED,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
})


def fig1_confound_scatter():
    df = pd.read_csv(f"{RESULTS_DIR}/axis_validation_per_cell_axis_position.csv")
    with open(f"{RESULTS_DIR}/axis_validation_summary.json") as f:
        summary = json.load(f)
    r = summary["pearson_r_axis_vs_n_genes"]

    fig, ax = plt.subplots(figsize=(6.0, 4.5))

    for age, color, label in [("Old", COLOR_OLD, "Old"), ("Young", COLOR_YOUNG, "Young")]:
        sub = df[df["age_group"] == age]
        ax.scatter(
            sub["n_genes_detected"], sub["axis_position"],
            s=14, color=color, alpha=0.55, linewidths=0, label=label,
        )

    ax.set_xlabel("Detected genes per cell")
    ax.set_ylabel("Position on Old→Young axis")
    ax.grid(True, color=COLOR_GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)

    ax.annotate(
        f"Pearson r = {r:.2f}",
        xy=(0.03, 0.95), xycoords="axes fraction",
        fontsize=11, color=COLOR_INK, va="top",
    )

    legend = ax.legend(frameon=False, loc="lower right", title=None)
    for lh in legend.legend_handles:
        lh.set_alpha(1.0)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = f"{FIG_DIR}/fig1_axis_vs_ngenes_confound.{ext}"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig2_deconfound_forest():
    with open(f"{RESULTS_DIR}/axis_validation_deconfound.json") as f:
        dec = json.load(f)

    rows = [
        ("Raw embeddings", dec["baseline"]),
        ("Detected-gene count\nregressed out", dec["approach1_residualized"]),
        ("Matched depth band", dec["approach2_matched_band"]),
    ]
    labels = [r[0] for r in rows]
    accuracy = [r[1]["loo_correct"] / r[1]["loo_n"] for r in rows]
    perm_p = [r[1]["perm_p"] for r in rows]
    y_pos = list(range(len(rows)))[::-1]  # top-to-bottom in listed order

    fig, (ax_acc, ax_p) = plt.subplots(1, 2, figsize=(9.0, 3.6), sharey=True)

    ax_acc.hlines(y_pos, 0, accuracy, color=COLOR_GRID, linewidth=1.5, zorder=1)
    ax_acc.scatter(accuracy, y_pos, s=90, color=COLOR_OLD, zorder=2)
    for y, v in zip(y_pos, accuracy):
        ax_acc.annotate(f"{v:.3f}", (v, y), xytext=(8, 0), textcoords="offset points",
                         va="center", fontsize=10, color=COLOR_INK)
    ax_acc.set_xlim(0, 1.05)
    ax_acc.set_xlabel("Donor LOO accuracy")
    ax_acc.set_yticks(y_pos, labels)
    ax_acc.grid(axis="x", color=COLOR_GRID, linewidth=0.6, zorder=0)
    ax_acc.set_axisbelow(True)

    ax_p.hlines(y_pos, 0, perm_p, color=COLOR_GRID, linewidth=1.5, zorder=1)
    ax_p.scatter(perm_p, y_pos, s=90, color=COLOR_OLD, zorder=2)
    for y, v in zip(y_pos, perm_p):
        ax_p.annotate(f"{v:.3f}", (v, y), xytext=(8, 0), textcoords="offset points",
                       va="center", fontsize=10, color=COLOR_INK)
    ax_p.set_xlim(0, max(perm_p) * 1.25)
    ax_p.set_ylim(min(y_pos) - 0.6, max(y_pos) + 0.8)
    ax_p.axvline(0.05, color=COLOR_REF, linewidth=1.2, linestyle="--", zorder=1)
    ax_p.annotate("p = 0.05", (0.05, max(y_pos) + 0.55), color=COLOR_REF,
                  fontsize=9, ha="center", va="bottom")
    ax_p.set_xlabel("Permutation p (donor-label shuffle)")
    ax_p.grid(axis="x", color=COLOR_GRID, linewidth=0.6, zorder=0)
    ax_p.set_axisbelow(True)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = f"{FIG_DIR}/fig2_deconfound_forest.{ext}"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig1_confound_scatter()
    fig2_deconfound_forest()
    print("Saved:")
    print(f"  {FIG_DIR}/fig1_axis_vs_ngenes_confound.pdf")
    print(f"  {FIG_DIR}/fig1_axis_vs_ngenes_confound.png")
    print(f"  {FIG_DIR}/fig2_deconfound_forest.pdf")
    print(f"  {FIG_DIR}/fig2_deconfound_forest.png")
