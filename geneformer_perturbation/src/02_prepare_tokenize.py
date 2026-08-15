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
IN_H5AD = "data/raw/facs/tms_facs_limb_msc_raw_counts.h5ad"
ORTHO_CSV = "data/mouse_to_human_orthologs.csv"
PREP_DIR = "data/prepped"              # where the cleaned .h5ad goes
TOK_DIR = "data/tokenized"             # tokenizer output
OUT_PREFIX = "tms_facs_limb_msc"
# columns in adata.obs to carry through for defining states and donor-level aggregation.
# Checked directly against the loaded FACS h5ad (PROGRESS.md 4.3b/4.3c):
#   - young/old label: "age_group" (values "Old" / "Young", capitalized)
#   - donor/mouse ID:  "mouse.id" (dot, NOT "mouse_id" like the Droplet cohort)
#   - sex:             "sex" (values "female" / "male")
# "condition" does not exist in this atlas cohort (it's not the reprogramming arm), so it has
# been dropped; the guard below silently skipped it either way.
KEEP_OBS = {
    "age_group": "age_group",
    "mouse.id": "mouse_id",   # renamed on output, dots are awkward as a downstream column/attr name
    "sex": "sex",
    "subtissue": "subtissue",
}  # {source_col: kept_name}

# FACS `subtissue` has 4 raw values (see PROGRESS.md 4.3c): 3 are naming variants of the same
# limb-muscle site, and 1 ("Muscle Diaphragm") is a different anatomical site found only in
# 2 old mice (18_46_F, 18_47_F) with no young counterpart, i.e. confounded with age group.
SUBTISSUE_LIMB_ALIASES = {
    "ForelimbandHindlimb": "limb",
    "Muscle forelimb and hindlimb": "limb",
    "forelimb and hindlimb": "limb",
}
SUBTISSUE_DIAPHRAGM_RAW = "Muscle Diaphragm"
SUBTISSUE_DIAPHRAGM_LABEL = "diaphragm"

# Primary analysis excludes diaphragm cells (confounded with age group, PROGRESS.md 4.3c
# decision). Set to False to run the diaphragm sensitivity check instead.
EXCLUDE_DIAPHRAGM = True
# ---------------------------------------------------------------------------


def prepare():
    Path(PREP_DIR).mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(IN_H5AD)

    # Geneformer wants RAW counts (it builds a rank encoding, not log-normalized).
    # If your object is already normalized, reload the raw layer instead.
    # e.g. adata = adata.raw.to_adata()

    # --- normalise subtissue labels and (by default) drop diaphragm cells ----
    if "subtissue" in adata.obs:
        adata.obs["subtissue"] = adata.obs["subtissue"].astype(str)
        adata.obs["subtissue"] = adata.obs["subtissue"].replace(SUBTISSUE_LIMB_ALIASES)
        adata.obs["subtissue"] = adata.obs["subtissue"].replace(
            {SUBTISSUE_DIAPHRAGM_RAW: SUBTISSUE_DIAPHRAGM_LABEL}
        )

        if EXCLUDE_DIAPHRAGM:
            n_before_filter = adata.n_obs
            adata = adata[adata.obs["subtissue"] != SUBTISSUE_DIAPHRAGM_LABEL].copy()
            print(
                f"excluded diaphragm cells: kept {adata.n_obs} / {n_before_filter} "
                "(EXCLUDE_DIAPHRAGM=True, see PROGRESS.md 4.3c)"
            )
        else:
            n_diaphragm = (adata.obs["subtissue"] == SUBTISSUE_DIAPHRAGM_LABEL).sum()
            print(
                f"EXCLUDE_DIAPHRAGM=False: keeping {n_diaphragm} diaphragm cells "
                "for the sensitivity run"
            )

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

    # model_version="V1" is required here: TranscriptomeTokenizer defaults to V2 settings
    # (gc104M token dictionary + gene median file, model_input_size 4096, special_token True),
    # which do NOT match our V1-10M checkpoint (gc30M dictionary, input size 2048, no special
    # tokens). Confirmed against examples/tokenizing_scRNAseq_data.ipynb cell 4, PROGRESS.md 4.4a.
    tk = TranscriptomeTokenizer(custom_attrs, nproc=4, model_version="V1")
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
