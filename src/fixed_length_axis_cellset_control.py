"""
Control for src/fixed_length_axis.py: 116 of the 118 cells dropped by the 1024-token cutoff
were Old (results/fixed_length_axis_summary.json drop_info), so the improved separation seen
there could reflect removing shallow Old cells rather than the length fix itself. This script
re-runs the ORIGINAL variable-length embedding extraction restricted to exactly the same 697
cells that survived the 1024-token cutoff (same cell set, no truncation), so the only
difference between this run and src/fixed_length_axis.py is the truncation step, not the cell
set, and the only difference between this run and the full-815-cell baseline
(results/axis_validation_summary.json) is the cell set, not the truncation.

Run:    python src/fixed_length_axis_cellset_control.py
Output: results/fixed_length_axis_cellset_control_cell_embeddings.csv
        results/fixed_length_axis_cellset_control_summary.json
Report only, no interpretation.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from datasets import load_from_disk
from scipy import sparse as sp

from fixed_length_axis import (
    ORIG_TOK_DATASET,
    PREPPED_H5AD,
    MODEL_DIR,
    OUT_DIR,
    FIXED_LENGTH,
    analyze,
)

VARLEN_TOK_DATASET = "data/tokenized/tms_facs_limb_msc_kept697_varlen.dataset"


def build_kept_variable_length_dataset():
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
    keep_mask = lengths >= FIXED_LENGTH   # SAME cell set as fixed_length_axis.py, no truncation
    n_kept = int(keep_mask.sum())
    kept_age_counts = pd.Series(age[keep_mask]).value_counts().to_dict()

    print(f"cells with tokenized length >= {FIXED_LENGTH} (same cell set as the truncated run): "
          f"{n_kept} / {len(ds)}")
    print(f"kept cells by age group: {kept_age_counts}")

    kept_ds = ds.select(np.where(keep_mask)[0].tolist())   # input_ids/length left untouched

    Path(VARLEN_TOK_DATASET).parent.mkdir(parents=True, exist_ok=True)
    kept_ds.save_to_disk(VARLEN_TOK_DATASET)
    print(f"wrote {VARLEN_TOK_DATASET} ({len(kept_ds)} cells, variable length, untruncated)")
    return n_kept


def extract_variable_length_embeddings(n_kept):
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
    embs_df = embex.extract_embs(
        MODEL_DIR,
        VARLEN_TOK_DATASET,
        OUT_DIR,
        "fixed_length_axis_cellset_control_raw_embs",
    )
    return embs_df


if __name__ == "__main__":
    Path(OUT_DIR).mkdir(exist_ok=True)

    n_kept = build_kept_variable_length_dataset()
    embs_df = extract_variable_length_embeddings(n_kept)

    embs_out = Path(OUT_DIR) / "fixed_length_axis_cellset_control_cell_embeddings.csv"
    embs_df.to_csv(embs_out, index=False)
    print(f"wrote {embs_out}")

    summary, donor_loo_df, axis_position = analyze(embs_df)
    summary["note"] = "same 697 cells as fixed_length_axis.py, but NO truncation (original variable length)"

    donor_loo_out = Path(OUT_DIR) / "fixed_length_axis_cellset_control_donor_loo.csv"
    donor_loo_df.to_csv(donor_loo_out, index=False)

    summary_out = Path(OUT_DIR) / "fixed_length_axis_cellset_control_summary.json"
    with open(summary_out, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"wrote {donor_loo_out}")
    print(f"wrote {summary_out}")
    print()
    print("=== donor LOO table ===")
    print(donor_loo_df.to_string(index=False))
    print()
    print("=== summary ===")
    print(json.dumps(summary, indent=2))
