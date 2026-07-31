"""
Token-length sweep for the Old-to-Young axis confound (follow-up to src/fixed_length_axis.py
and src/fixed_length_axis_cellset_control.py). Runs both the fixed-length-truncation condition
and the same-cell-set variable-length control condition across a range of token-length cutoffs,
to see how donor separation and the detected-gene-count correlation move as the cutoff moves,
and to what extent that movement is driven by truncation itself versus which cells the cutoff
drops.

At each length L in LENGTHS:
  fixed:            keep cells with tokenized length >= L, truncate their input_ids to the
                     first L (highest-ranked) tokens, set length=L.
  cellset_control:   keep the SAME cells (length >= L), but leave their input_ids/length
                     untouched (original variable length). Isolates the cell-drop effect from
                     the truncation effect, same logic as fixed_length_axis_cellset_control.py.

At L=2048 (the V1-10M input cap, and the max length any cell reaches) both conditions keep
only cells with length exactly 2048, and truncating a length-2048 cell to 2048 tokens is a
no-op, so fixed and cellset_control necessarily converge to identical numbers there. That is
expected and is not special-cased away, it is a sanity check on the pipeline.

Everything else matches the prior runs exactly: ortholog-mapped, diaphragm-excluded FACS data
(data/prepped/tms_facs_limb_msc.h5ad), Geneformer V1-10M, same cell-level global-centroid /
donor-mean leave-one-donor-out procedure, 1000-draw donor-label permutation with seed=42, and
pre-truncation (original) detected-gene count as the covariate for the Pearson correlation.

Run:    python src/token_length_sweep.py
Output: results/token_length_sweep.csv   (tidy long format, one row per length x condition)
Report only, no plotting, no interpretation.
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
SWEEP_TOK_DIR = "data/tokenized/sweep"

MODEL_DIR = "Geneformer/Geneformer-V1-10M"
OUT_DIR = "results"

LENGTHS = [256, 384, 512, 768, 1024, 1536, 2048]
STATE_KEY = "age_group"
START_STATE = "Old"
GOAL_STATE = "Young"

SEED = 42
N_PERM = 1000
# ---------------------------------------------------------------------------


def load_base_dataset_with_ngenes():
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
    return ds


def cos_sim(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def donor_loo_accuracy(donor_means, donor_age):
    donors = list(donor_means.keys())
    n = len(donors)
    correct = 0
    for held_out in donors:
        others = [d for d in donors if d != held_out]
        old_others = [d for d in others if donor_age[d] == START_STATE]
        young_others = [d for d in others if donor_age[d] == GOAL_STATE]
        o_cent = np.mean([donor_means[d] for d in old_others], axis=0)
        y_cent = np.mean([donor_means[d] for d in young_others], axis=0)
        d_old = cos_sim(donor_means[held_out], o_cent)
        d_young = cos_sim(donor_means[held_out], y_cent)
        pred = GOAL_STATE if d_young > d_old else START_STATE
        correct += int(pred == donor_age[held_out])
    return correct, n


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

    old_centroid = E[old_mask].mean(axis=0)
    young_centroid = E[young_mask].mean(axis=0)
    centroid_dist = float(1 - cos_sim(old_centroid, young_centroid))

    axis_position = np.array([
        cos_sim(E[i], young_centroid) - cos_sim(E[i], old_centroid) for i in range(n_cells)
    ])
    r, p = pearsonr(axis_position, ngenes)

    donor_means = {}
    donor_age = {}
    for m in np.unique(mouse):
        mask = mouse == m
        donor_means[m] = E[mask].mean(axis=0)
        donor_age[m] = age[mask][0]

    loo_correct, loo_n = donor_loo_accuracy(donor_means, donor_age)

    rng = np.random.default_rng(SEED)
    donors = list(donor_means.keys())
    true_ages = np.array([donor_age[d] for d in donors])
    perm_correct = np.empty(N_PERM, dtype=int)
    for i in range(N_PERM):
        shuffled_ages = rng.permutation(true_ages)
        shuffled_age_map = dict(zip(donors, shuffled_ages))
        p_correct, _ = donor_loo_accuracy(donor_means, shuffled_age_map)
        perm_correct[i] = p_correct
    perm_p = float((perm_correct >= loo_correct).sum()) / N_PERM

    return {
        "loo_accuracy": loo_correct / loo_n,
        "loo_correct": loo_correct,
        "loo_n": loo_n,
        "perm_p": perm_p,
        "pearson_r": float(r),
        "pearson_p": float(p),
        "centroid_dist": centroid_dist,
    }


def run_condition(base_ds, length, condition, prepped_ages):
    lengths = np.array(base_ds["length"])
    age = np.array(base_ds["age_group"])
    keep_mask = lengths >= length
    n_kept = int(keep_mask.sum())
    dropped_age_counts = pd.Series(age[~keep_mask]).value_counts().to_dict()
    n_dropped_old = int(dropped_age_counts.get(START_STATE, 0))
    n_dropped_young = int(dropped_age_counts.get(GOAL_STATE, 0))

    kept_ds = base_ds.select(np.where(keep_mask)[0].tolist())

    if condition == "fixed":
        def truncate(example):
            example["input_ids"] = example["input_ids"][:length]
            example["length"] = length
            return example
        kept_ds = kept_ds.map(truncate)

    tok_path = f"{SWEEP_TOK_DIR}/{condition}_L{length}.dataset"
    Path(tok_path).parent.mkdir(parents=True, exist_ok=True)
    kept_ds.save_to_disk(tok_path)

    from geneformer import EmbExtractor
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
    embs_df = embex.extract_embs(MODEL_DIR, tok_path, OUT_DIR, f"sweep_{condition}_L{length}_embs")

    stats = analyze(embs_df)

    print(f"L={length} condition={condition}: kept={n_kept} "
          f"dropped_old={n_dropped_old} dropped_young={n_dropped_young} "
          f"loo={stats['loo_correct']}/{stats['loo_n']} perm_p={stats['perm_p']:.4f} "
          f"r={stats['pearson_r']:.4f} centroid_dist={stats['centroid_dist']:.6f}")

    return {
        "token_length": length,
        "condition": condition,
        "n_kept": n_kept,
        "n_dropped_old": n_dropped_old,
        "n_dropped_young": n_dropped_young,
        "loo_accuracy": stats["loo_accuracy"],
        "loo_correct": stats["loo_correct"],
        "loo_n": stats["loo_n"],
        "perm_p": stats["perm_p"],
        "pearson_r": stats["pearson_r"],
        "pearson_p": stats["pearson_p"],
        "centroid_dist": stats["centroid_dist"],
    }


if __name__ == "__main__":
    Path(OUT_DIR).mkdir(exist_ok=True)

    base_ds = load_base_dataset_with_ngenes()
    prepped_ages = None  # unused, kept for signature symmetry

    rows = []
    for length in LENGTHS:
        for condition in ["fixed", "cellset_control"]:
            row = run_condition(base_ds, length, condition, prepped_ages)
            rows.append(row)

    sweep_df = pd.DataFrame(rows)
    tidy_cols = ["token_length", "condition", "n_kept", "n_dropped_old", "n_dropped_young",
                 "loo_accuracy", "perm_p", "pearson_r", "centroid_dist"]
    out_path = Path(OUT_DIR) / "token_length_sweep.csv"
    sweep_df[tidy_cols].to_csv(out_path, index=False)

    full_out_path = Path(OUT_DIR) / "token_length_sweep_full.csv"
    sweep_df.to_csv(full_out_path, index=False)

    print(f"\nwrote {out_path}")
    print(f"wrote {full_out_path} (includes loo_correct/loo_n/pearson_p for reference)")
    print()
    print(sweep_df[tidy_cols].to_string(index=False))
