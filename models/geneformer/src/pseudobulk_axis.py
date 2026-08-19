"""
Experiment (supervisor feedback): build the Old-to-Young axis from donor-level pseudobulk
rather than individual cells, and compare donor separation against the cell-level baseline.

Steps:
  1. Aggregate each donor's raw counts (from data/prepped/tms_facs_limb_msc.h5ad, the same
     ortholog-mapped, diaphragm-excluded prepped file the cell-level pipeline uses) into one
     pseudobulk profile per donor by summing raw counts across that donor's cells. 14 donors
     (8 Old, 6 Young), same as the cell-level donor-level analyses elsewhere in this project.
  2. Tokenize the 14 pseudobulk profiles with the identical V1 tokenizer settings as
     src/02_prepare_tokenize.py.
  3. Extract a per-donor embedding for each of the 14 profiles with Geneformer (same V1-10M
     checkpoint and settings as src/03_in_silico_perturb.py's get_state_embs, but per-profile
     rather than a state mean: summary_stat=None).
  4. Form Old/Young centroids from the 14 donor embeddings, compute each donor's axis_position
     (cos_sim to Young centroid minus cos_sim to Old centroid, same convention used at cell
     level), and run leave-one-donor-out (LODO) nearest-centroid classification: for each
     donor, build Old/Young centroids from the OTHER 13 donors only, classify the held-out
     donor by nearest centroid, tally correct.
  5. Correlate donor axis_position against donor-level detected-gene count. Detected-gene
     count here is the SAME quantity already established at donor level elsewhere in this
     project (results/axis_validation_donor_ngenes_distribution.csv): the per-donor mean of
     each cell's number of nonzero-count genes, NOT the number of nonzero genes in the
     pseudobulk profile itself (which would be close to all genes for every donor after
     summing dozens of cells, and would not be a comparable quantity to the cell-level
     n_genes_detected this is being checked against).

Run:    python src/pseudobulk_axis.py
Output: results/pseudobulk_axis_donor_embeddings.csv
        results/pseudobulk_axis_summary.json
Report only, no interpretation (this script does not draw conclusions, callers should not
either until asked).
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse as sp
from scipy.stats import pearsonr

# ---------------------------------------------------------------- TODO -----
PREPPED_H5AD = "data/prepped/tms_facs_limb_msc.h5ad"          # cell-level prepped file
# own directory, not data/prepped/: TranscriptomeTokenizer.tokenize_data scans its entire
# data_directory for every .h5ad file it contains (same directory-wide-scan gotcha as
# InSilicoPerturberStats, see 03_in_silico_perturb.py module docstring) -- sharing
# data/prepped/ with tms_facs_limb_msc.h5ad would silently tokenize both together.
PSEUDOBULK_PREP_DIR = "data/prepped_pseudobulk"
PSEUDOBULK_H5AD = f"{PSEUDOBULK_PREP_DIR}/tms_facs_limb_msc_pseudobulk.h5ad"
TOK_DIR = "data/tokenized"
TOK_PREFIX = "tms_facs_limb_msc_pseudobulk"
TOK_DATASET = f"{TOK_DIR}/{TOK_PREFIX}.dataset"

MODEL_DIR = "Geneformer/Geneformer-V1-10M"

OUT_DIR = "results"
DONOR_NGENES_CSV = "results/axis_validation_donor_ngenes_distribution.csv"  # already computed

STATE_KEY = "age_group"
START_STATE = "Old"
GOAL_STATE = "Young"
# ---------------------------------------------------------------------------


def build_pseudobulk():
    adata = sc.read_h5ad(PREPPED_H5AD)
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X)

    donors = adata.obs["mouse_id"].astype(str)
    donor_order = sorted(donors.unique())

    pb_rows = []
    obs_rows = []
    for donor in donor_order:
        mask = (donors == donor).values
        pb_rows.append(X[mask].sum(axis=0))

        sub = adata.obs.loc[mask]
        age_vals = sub["age_group"].unique()
        sex_vals = sub["sex"].unique()
        subtissue_vals = sub["subtissue"].unique()
        assert len(age_vals) == 1, f"donor {donor} has mixed age_group: {age_vals}"
        assert len(sex_vals) == 1, f"donor {donor} has mixed sex: {sex_vals}"
        obs_rows.append({
            "mouse_id": donor,
            "age_group": age_vals[0],
            "sex": sex_vals[0],
            "subtissue": subtissue_vals[0],
            "n_cells_aggregated": int(mask.sum()),
        })

    pb_X = np.vstack(pb_rows).astype(np.float32)
    obs_df = pd.DataFrame(obs_rows).set_index("mouse_id", drop=False)

    pb_adata = sc.AnnData(X=pb_X, obs=obs_df, var=adata.var.copy())
    pb_adata.obs["n_counts"] = pb_adata.X.sum(axis=1)

    print(f"pseudobulk profiles: {pb_adata.n_obs} donors x {pb_adata.n_vars} genes")
    print(obs_df[["age_group", "sex", "subtissue", "n_cells_aggregated"]].to_string())

    Path(PSEUDOBULK_H5AD).parent.mkdir(parents=True, exist_ok=True)
    pb_adata.write_h5ad(PSEUDOBULK_H5AD)
    print(f"wrote {PSEUDOBULK_H5AD}")
    return PSEUDOBULK_H5AD


def tokenize_pseudobulk():
    from geneformer import TranscriptomeTokenizer

    custom_attrs = {"mouse_id": "mouse_id", "age_group": "age_group"}
    # same V1 settings as src/02_prepare_tokenize.py (PROGRESS.md 4.4a): TranscriptomeTokenizer
    # defaults to V2/gc104M otherwise.
    tk = TranscriptomeTokenizer(custom_attrs, nproc=1, model_version="V1")
    tk.tokenize_data(
        data_directory=PSEUDOBULK_PREP_DIR,
        output_directory=TOK_DIR,
        output_prefix=TOK_PREFIX,
        file_format="h5ad",
    )
    print(f"tokenized pseudobulk dataset written to {TOK_DATASET}")


def extract_donor_embeddings():
    from geneformer import EmbExtractor

    embex = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cell",
        emb_label=["mouse_id", "age_group"],
        max_ncells=14,             # all 14 donor profiles, no subsampling
        emb_layer=0,
        summary_stat=None,         # per-profile embeddings, not a state mean
        forward_batch_size=14,
        model_version="V1",        # gc30M dictionary, matches V1-10M
        nproc=1,
    )
    embs_df = embex.extract_embs(
        MODEL_DIR,
        TOK_DATASET,
        OUT_DIR,
        "pseudobulk_axis_raw_embs",
    )
    return embs_df


def cos_sim(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def analyze(embs_df):
    emb_cols = [c for c in embs_df.columns if c not in ("mouse_id", "age_group")]
    embs_df = embs_df.reset_index(drop=True)
    E = embs_df[emb_cols].to_numpy(dtype=np.float64)
    age = embs_df["age_group"].to_numpy()
    donors = embs_df["mouse_id"].to_numpy()
    n = len(embs_df)

    old_mask = age == START_STATE
    young_mask = age == GOAL_STATE
    assert old_mask.sum() + young_mask.sum() == n, "unexpected age_group values"

    # --- global centroids (all 14 donors) and axis_position, same convention as cell level:
    #     axis_position = cos_sim(donor, Young_centroid) - cos_sim(donor, Old_centroid)
    old_centroid = E[old_mask].mean(axis=0)
    young_centroid = E[young_mask].mean(axis=0)
    # 1 - cos_sim, same "distance" convention as old_spread/young_spread below and as
    # results/axis_validation_summary.json's centroid_cos_dist at cell level
    centroid_cos_dist = float(1 - cos_sim(old_centroid, young_centroid))

    axis_position = np.array([
        cos_sim(E[i], young_centroid) - cos_sim(E[i], old_centroid) for i in range(n)
    ])

    old_spread = np.array([cos_sim(E[i], old_centroid) for i in range(n) if old_mask[i]])
    young_spread = np.array([cos_sim(E[i], young_centroid) for i in range(n) if young_mask[i]])
    # spread = 1 - mean cos_sim to own centroid (dispersion), same convention as cell-level report
    old_spread_mean = float(1 - old_spread.mean())
    young_spread_mean = float(1 - young_spread.mean())

    # --- leave-one-donor-out nearest-centroid classification -----------------------------
    loo_correct = 0
    loo_rows = []
    for i in range(n):
        others = np.arange(n) != i
        o_cent = E[others & old_mask].mean(axis=0)
        y_cent = E[others & young_mask].mean(axis=0)
        d_old = cos_sim(E[i], o_cent)
        d_young = cos_sim(E[i], y_cent)
        pred = GOAL_STATE if d_young > d_old else START_STATE
        correct = pred == age[i]
        loo_correct += int(correct)
        loo_rows.append({
            "mouse_id": donors[i], "age_group": age[i], "predicted": pred,
            "correct": correct, "cos_to_old_centroid": d_old, "cos_to_young_centroid": d_young,
        })

    loo_df = pd.DataFrame(loo_rows)

    # --- correlation with donor-level detected-gene count --------------------------------
    ngenes_df = pd.read_csv(DONOR_NGENES_CSV)[["mouse_id", "mean"]].rename(
        columns={"mean": "mean_n_genes_detected"}
    )
    merged = pd.DataFrame({"mouse_id": donors, "axis_position": axis_position}).merge(
        ngenes_df, on="mouse_id", how="left"
    )
    assert merged["mean_n_genes_detected"].notna().all(), "donor ID mismatch with ngenes CSV"

    r, p = pearsonr(merged["axis_position"], merged["mean_n_genes_detected"])

    summary = {
        "n_donors": n,
        "n_old_donors": int(old_mask.sum()),
        "n_young_donors": int(young_mask.sum()),
        "centroid_cos_dist": centroid_cos_dist,
        "old_spread_mean": old_spread_mean,
        "young_spread_mean": young_spread_mean,
        "loo_correct": loo_correct,
        "loo_n": n,
        "pearson_r_axis_vs_donor_n_genes": float(r),
        "pearson_p_axis_vs_donor_n_genes": float(p),
        "cell_level_baseline_loo": "11/14",
        "cell_level_baseline_pearson_r": 0.6909679352631691,
    }

    age_lookup = dict(zip(donors, age))
    donor_table = merged.merge(loo_df[["mouse_id", "predicted", "correct"]], on="mouse_id")
    donor_table["age_group"] = donor_table["mouse_id"].map(age_lookup)
    donor_table = donor_table[["mouse_id", "age_group", "axis_position", "mean_n_genes_detected", "predicted", "correct"]]
    donor_table = donor_table.sort_values(["age_group", "mouse_id"]).reset_index(drop=True)

    return summary, donor_table


if __name__ == "__main__":
    Path(OUT_DIR).mkdir(exist_ok=True)

    build_pseudobulk()
    tokenize_pseudobulk()
    embs_df = extract_donor_embeddings()

    embs_out = Path(OUT_DIR) / "pseudobulk_axis_donor_embeddings.csv"
    embs_df.to_csv(embs_out, index=False)
    print(f"wrote {embs_out}")

    summary, donor_table = analyze(embs_df)

    donor_table_out = Path(OUT_DIR) / "pseudobulk_axis_donor_table.csv"
    donor_table.to_csv(donor_table_out, index=False)

    summary_out = Path(OUT_DIR) / "pseudobulk_axis_summary.json"
    with open(summary_out, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"wrote {donor_table_out}")
    print(f"wrote {summary_out}")
    print()
    print("=== per-donor table ===")
    print(donor_table.to_string(index=False))
    print()
    print("=== summary ===")
    print(json.dumps(summary, indent=2))
