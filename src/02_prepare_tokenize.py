"""
Step 2. Prepare AnnData and tokenize for Geneformer.

Takes a mouse scRNA-seq .h5ad (limb-muscle MSC subset of TMS, or GSE176206),
maps genes to human Ensembl IDs, adds the two fields the Geneformer tokenizer
needs (`ensembl_id` in var, `n_counts` in obs), keeps the metadata columns we'll
need later (age group, condition, timepoint), and runs the tokenizer.

Run:     python src/02_prepare_tokenize.py
Output:  data/tokenized/  (a Hugging Face dataset the perturber reads)

NOTE: confirm the tokenizer arguments against the example notebook in the
Geneformer version you cloned (examples/tokenizing_scRNAseq_data.ipynb).
Arg names have changed between V1 and V2.
"""

from pathlib import Path
import scanpy as sc
import pandas as pd

# ---------------------------------------------------------------- TODO -----
IN_H5AD = "data/raw/tms_droplet_limb_msc_raw_counts.h5ad"
ORTHO_CSV = "data/mouse_to_human_orthologs.csv"
PREP_DIR = "data/prepped"              # where the cleaned .h5ad goes
TOK_DIR = "data/tokenized"             # tokenizer output
OUT_PREFIX = "tms_droplet_limb_msc"
# columns in adata.obs to carry through so we can define old/young + condition later
KEEP_OBS = {"age": "age", "condition": "condition"}  # {source_col: kept_name}
# ---------------------------------------------------------------------------


def prepare():
    Path(PREP_DIR).mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(IN_H5AD)

    # Geneformer wants RAW counts (it builds a rank encoding, not log-normalized).
    # If your object is already normalized, reload the raw layer instead.
    # e.g. adata = adata.raw.to_adata()

    # --- map mouse genes -> human Ensembl ------------------------------------
    ortho = pd.read_csv(ORTHO_CSV)
    sym2ens = dict(zip(ortho["mouse_symbol"], ortho["human_ensembl"]))

    # assumes var_names are mouse gene symbols; adjust if yours are mouse Ensembl
    adata.var["human_ensembl"] = adata.var_names.map(sym2ens)
    n_before = adata.n_vars
    adata = adata[:, adata.var["human_ensembl"].notna()].copy()
    print(f"kept {adata.n_vars} / {n_before} genes with a human ortholog")

    # tokenizer reads gene ids from var["ensembl_id"]
    adata.var["ensembl_id"] = adata.var["human_ensembl"].values

    # tokenizer needs total counts per cell
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)

    # keep the metadata we'll need for defining states
    for src, keep in KEEP_OBS.items():
        if src in adata.obs:
            adata.obs[keep] = adata.obs[src].astype(str)

    out = Path(PREP_DIR) / f"{OUT_PREFIX}.h5ad"
    adata.write_h5ad(out)
    print(f"wrote {out}")
    return PREP_DIR


def tokenize(prep_dir):
    from geneformer import TranscriptomeTokenizer

    # keys = obs columns to keep in the tokenized dataset (value = name in output)
    custom_attrs = {name: name for name in KEEP_OBS.values()}

    tk = TranscriptomeTokenizer(custom_attrs, nproc=4)
    tk.tokenize_data(
        data_directory=prep_dir,   # dir containing the .h5ad
        output_directory=TOK_DIR,
        output_prefix=OUT_PREFIX,
        file_format="h5ad",
    )
    print(f"tokenized dataset written under {TOK_DIR}/")


if __name__ == "__main__":
    prep_dir = prepare()
    tokenize(prep_dir)
