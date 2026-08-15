"""
Experiment: re-extract per-cell embeddings after truncating every cell's tokenized input to a
fixed 1024 top-ranked genes, discarding cells whose tokenized length is below 1024, then rebuild
the Old/Young axis on these fixed-length embeddings and compare against the original (untruncated,
variable-length) cell-level results.

Why: the original per-cell embeddings correlate with detected-gene count (n_genes_detected) at
r=0.69 (results/axis_validation_summary.json). Geneformer's rank-value encoding gives each cell a
token sequence as long as its own number of detected genes (capped at 2048 for V1). Forcing every
retained cell to the SAME fixed input length (1024 top-ranked genes) removes sequence-length as a
source of cell-to-cell difference in what the model sees, as a further, more direct deconfounding
check alongside the two already done (embedding residualization, matched detected-gene-count band;
PROGRESS.md 6.6).

Alignment note: the on-disk tokenized dataset (data/tokenized/tms_facs_limb_msc.dataset) preserves
the exact row order of data/prepped/tms_facs_limb_msc.h5ad (verified positionally, all 815 rows
identical mouse_id/age_group order) -- TranscriptomeTokenizer does not shuffle at tokenize time,
only EmbExtractor's internal downsample_and_sort does, at extraction time. n_genes_detected is
therefore computed once here directly from the prepped h5ad's raw counts (nonzero genes on the
same 14314-gene ortholog-mapped set used for tokenization) and attached to the dataset BEFORE
extraction, as a custom label column alongside mouse_id/age_group, so it survives extract_embs's
internal reordering correctly instead of relying on fragile positional matching against a
separately-generated CSV.

Run:    python src/fixed_length_axis.py
Output: results/fixed_length_axis_cell_embeddings.csv
        results/fixed_length_axis_summary.json
Report only, no interpretation.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from datasets import load_from_disk
from scipy import sparse as sp
from scipy.stats import pearsonr

# ---------------------------------------------------------------- TODO -----
PREPPED_H5AD = "data/prepped/tms_facs_limb_msc.h5ad"
ORIG_TOK_DATASET = "data/tokenized/tms_facs_limb_msc.dataset"
TRUNC_TOK_DATASET = "data/tokenized/tms_facs_limb_msc_trunc1024.dataset"

MODEL_DIR = "Geneformer/Geneformer-V1-10M"
OUT_DIR = "results"

FIXED_LENGTH = 1024
STATE_KEY = "age_group"
START_STATE = "Old"
GOAL_STATE = "Young"

SEED = 42
N_PERM = 1000
# ---------------------------------------------------------------------------


def build_truncated_dataset():
    ds = load_from_disk(ORIG_TOK_DATASET)

    prepped = sc.read_h5ad(PREPPED_H5AD)
    assert prepped.n_obs == len(ds), "row count mismatch between prepped h5ad and tokenized dataset"
    assert list(prepped.obs["mouse_id"].astype(str)) == list(ds["mouse_id"]), \
        "prepped h5ad and tokenized dataset are not positionally aligned"

    X = prepped.X
    if sp.issparse(X):
        X = X.toarray()
    n_genes_detected = (np.asarray(X) > 0).sum(axis=1).astype(int)
    ds = ds.add_column("n_genes_detected_orig", n_genes_detected.tolist())

    lengths = np.array(ds["length"])
    age = np.array(ds["age_group"])
    drop_mask = lengths < FIXED_LENGTH
    n_dropped = int(drop_mask.sum())
    dropped_age_counts = pd.Series(age[drop_mask]).value_counts().to_dict()

    print(f"cells with tokenized length < {FIXED_LENGTH}: {n_dropped} / {len(ds)}")
    print(f"dropped cells by age group: {dropped_age_counts}")

    kept_ds = ds.select(np.where(~drop_mask)[0].tolist())

    def truncate(example):
        example["input_ids"] = example["input_ids"][:FIXED_LENGTH]
        example["length"] = FIXED_LENGTH
        return example

    kept_ds = kept_ds.map(truncate)

    Path(TRUNC_TOK_DATASET).parent.mkdir(parents=True, exist_ok=True)
    kept_ds.save_to_disk(TRUNC_TOK_DATASET)
    print(f"wrote {TRUNC_TOK_DATASET} ({len(kept_ds)} cells)")

    drop_info = {
        "n_total_cells": len(ds),
        "n_dropped": n_dropped,
        "n_kept": len(kept_ds),
        "dropped_age_counts": dropped_age_counts,
        "fixed_length": FIXED_LENGTH,
    }
    return drop_info


def extract_truncated_embeddings():
    from geneformer import EmbExtractor

    kept_ds = load_from_disk(TRUNC_TOK_DATASET)
    n_kept = len(kept_ds)

    embex = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cell",
        emb_label=["mouse_id", "age_group", "n_genes_detected_orig"],
        max_ncells=n_kept,
        emb_layer=0,
        summary_stat=None,
        forward_batch_size=16,
        model_version="V1",
        nproc=4,
    )
    embs_df = embex.extract_embs(
        MODEL_DIR,
        TRUNC_TOK_DATASET,
        OUT_DIR,
        "fixed_length_axis_raw_embs",
    )
    return embs_df


def cos_sim(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def donor_loo_accuracy(donor_means, donor_age):
    """donor_means: dict mouse_id -> embedding vector. donor_age: dict mouse_id -> age label.
    Returns (n_correct, n_donors, per-donor prediction rows)."""
    donors = list(donor_means.keys())
    n = len(donors)
    correct = 0
    rows = []
    for i, held_out in enumerate(donors):
        others = [d for d in donors if d != held_out]
        old_others = [d for d in others if donor_age[d] == START_STATE]
        young_others = [d for d in others if donor_age[d] == GOAL_STATE]
        o_cent = np.mean([donor_means[d] for d in old_others], axis=0)
        y_cent = np.mean([donor_means[d] for d in young_others], axis=0)
        d_old = cos_sim(donor_means[held_out], o_cent)
        d_young = cos_sim(donor_means[held_out], y_cent)
        pred = GOAL_STATE if d_young > d_old else START_STATE
        is_correct = pred == donor_age[held_out]
        correct += int(is_correct)
        rows.append({"mouse_id": held_out, "age_group": donor_age[held_out],
                      "predicted": pred, "correct": is_correct})
    return correct, n, rows


def analyze(embs_df):
    emb_cols = [c for c in embs_df.columns if c not in ("mouse_id", "age_group", "n_genes_detected_orig")]
    embs_df = embs_df.reset_index(drop=True)
    E = embs_df[emb_cols].to_numpy(dtype=np.float64)
    age = embs_df["age_group"].to_numpy()
    mouse = embs_df["mouse_id"].to_numpy()
    ngenes = embs_df["n_genes_detected_orig"].to_numpy(dtype=float)
    n_cells = len(embs_df)

    old_mask = age == START_STATE
    young_mask = age == GOAL_STATE

    # --- cell-level global centroids (mirrors results/axis_validation_summary.json) --------
    old_centroid = E[old_mask].mean(axis=0)
    young_centroid = E[young_mask].mean(axis=0)
    centroid_cos_dist = float(1 - cos_sim(old_centroid, young_centroid))

    axis_position = np.array([
        cos_sim(E[i], young_centroid) - cos_sim(E[i], old_centroid) for i in range(n_cells)
    ])
    old_spread = np.array([cos_sim(E[i], old_centroid) for i in range(n_cells) if old_mask[i]])
    young_spread = np.array([cos_sim(E[i], young_centroid) for i in range(n_cells) if young_mask[i]])
    old_spread_mean = float(1 - old_spread.mean())
    young_spread_mean = float(1 - young_spread.mean())

    r, p = pearsonr(axis_position, ngenes)

    # --- donor-level leave-one-out nearest-centroid classification, on donor MEAN cell
    #     embeddings (kept cells only per donor), matching the original donor-level method ---
    donor_means = {}
    donor_age = {}
    for m in np.unique(mouse):
        mask = mouse == m
        donor_means[m] = E[mask].mean(axis=0)
        donor_age[m] = age[mask][0]

    loo_correct, loo_n, loo_rows = donor_loo_accuracy(donor_means, donor_age)

    # --- donor-label permutation control (1000 shuffles, seed=42), same procedure as
    #     results/axis_validation_permutation_accuracies.csv -----------------------------
    rng = np.random.default_rng(SEED)
    donors = list(donor_means.keys())
    true_ages = np.array([donor_age[d] for d in donors])
    perm_accuracies = []
    for _ in range(N_PERM):
        shuffled_ages = rng.permutation(true_ages)
        shuffled_age_map = dict(zip(donors, shuffled_ages))
        p_correct, _, _ = donor_loo_accuracy(donor_means, shuffled_age_map)
        perm_accuracies.append(p_correct)
    perm_accuracies = np.array(perm_accuracies)
    perm_ge_true_count = int((perm_accuracies >= loo_correct).sum())
    perm_p = perm_ge_true_count / N_PERM

    summary = {
        "fixed_length": FIXED_LENGTH,
        "n_cells": n_cells,
        "n_old_cells": int(old_mask.sum()),
        "n_young_cells": int(young_mask.sum()),
        "n_donors": len(donors),
        "centroid_cos_dist": centroid_cos_dist,
        "old_spread_mean": old_spread_mean,
        "young_spread_mean": young_spread_mean,
        "loo_correct": loo_correct,
        "loo_n": loo_n,
        "perm_mean": float(perm_accuracies.mean()),
        "perm_sd": float(perm_accuracies.std(ddof=1)),
        "perm_ge_true_count": perm_ge_true_count,
        "n_perm": N_PERM,
        "perm_seed": SEED,
        "perm_p": perm_p,
        "pearson_r_axis_vs_n_genes_detected_orig": float(r),
        "pearson_p_axis_vs_n_genes_detected_orig": float(p),
    }

    donor_loo_df = pd.DataFrame(loo_rows).sort_values(["age_group", "mouse_id"]).reset_index(drop=True)

    return summary, donor_loo_df, axis_position


if __name__ == "__main__":
    Path(OUT_DIR).mkdir(exist_ok=True)

    drop_info = build_truncated_dataset()
    embs_df = extract_truncated_embeddings()

    embs_out = Path(OUT_DIR) / "fixed_length_axis_cell_embeddings.csv"
    embs_df.to_csv(embs_out, index=False)
    print(f"wrote {embs_out}")

    summary, donor_loo_df, axis_position = analyze(embs_df)
    summary["drop_info"] = drop_info

    donor_loo_out = Path(OUT_DIR) / "fixed_length_axis_donor_loo.csv"
    donor_loo_df.to_csv(donor_loo_out, index=False)

    summary_out = Path(OUT_DIR) / "fixed_length_axis_summary.json"
    with open(summary_out, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"wrote {donor_loo_out}")
    print(f"wrote {summary_out}")
    print()
    print("=== drop info ===")
    print(json.dumps(drop_info, indent=2))
    print()
    print("=== donor LOO table ===")
    print(donor_loo_df.to_string(index=False))
    print()
    print("=== summary ===")
    print(json.dumps(summary, indent=2))
