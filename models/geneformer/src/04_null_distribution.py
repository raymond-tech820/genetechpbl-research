"""
Step 4. Random-gene null distribution for the perturbation ranking.

Runs the identical single-gene overexpression pipeline used for the 7 candidate factors
(src/03_in_silico_perturb.py, MODE="individual") on N randomly chosen genes instead, so the
candidates' goal-state shifts can eventually be judged against "what does a random gene do"
rather than read at face value.

Run:     python src/04_null_distribution.py
Output:  results/tms_facs_limb_msc_null_n<N>/  (per-gene raw pickles + a combined shifts CSV)

WHY EXPRESSION-MATCHED, NOT UNIFORM RANDOM:
Geneformer's rank-value encoding ranks each gene's expression against every other gene
detected in that cell. Overexpressing a gene moves it to the front of that ranking. A gene
that was barely detected has much further to move than a gene that was already highly
expressed, for purely mechanical/positional reasons unrelated to any biological effect. A
uniform-random null would therefore be comparing candidates (which may sit at any expression
level) against a null sample that likely sits at a different, uncontrolled expression level,
which is not a fair comparison. The fix mirrors two things already used elsewhere on this
project (PROGRESS.md 4.2, Kei's risk score background correction, and Zihan's expression-matched
null): bin all genes by expression level and only compare within the same bin.

This reuses `get_state_embs` and `run_one_perturbation` directly from 03_in_silico_perturb.py
(loaded by file path since its filename starts with a digit and can't be `import`ed normally)
rather than re-implementing them, so this null run is guaranteed to use IDENTICAL settings
(model, state definition, cell counts, model_version, filter_data, etc.) to the candidate runs,
not just similar-looking settings that could quietly drift apart.
"""

import importlib.util
import pickle
import random
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------- TODO -----
SCRIPT_03_PATH = Path(__file__).parent / "03_in_silico_perturb.py"

# Same prepped h5ad that 03's tokenized dataset was built from (post ortholog-mapping and
# diaphragm exclusion, PROGRESS.md 4.3c/4.4): this is where "expression level in our FACS
# data" is measured from, NOT Geneformer's pretrained Genecorpus gene-median file.
PREPPED_H5AD = "data/prepped/tms_facs_limb_msc.h5ad"

# Same V1-10M token dictionary used everywhere else in this pipeline (PROGRESS.md 3, 4.4a).
TOKEN_DICT_PATH = "Geneformer/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"

# OCT4 and NANOG were excluded from data/mouse_to_human_orthologs.csv entirely because they're
# one2many orthologs (PROGRESS.md 4.4a), NOT because they're unexpressed, so they're absent
# from PREPPED_H5AD and from the tokenized vocabulary. Their true mouse-gene expression is
# looked up directly here (by mouse symbol, from the raw FACS h5ad) purely so they still get
# placed in the correct expression bin; they are never draw-eligible as null genes themselves,
# only the unambiguous 1:1-ortholog genes we can actually perturb cleanly are.
RAW_FACS_H5AD = "data/raw/facs/tms_facs_limb_msc_raw_counts.h5ad"
MANUAL_MOUSE_SYMBOL = {
    "ENSG00000204531": "Pou5f1",  # OCT4
    "ENSG00000111704": "Nanog",   # NANOG
}

OUT_DIR = "results"
OUT_PREFIX = "tms_facs_limb_msc_null"

N_RANDOM_GENES = 50   # timing/correctness check size; scale up only after this is reviewed
N_EXPRESSION_BINS = 25   # matches Kei's sc.tl.score_genes n_bins=25 (PROGRESS.md 4.2) exactly;
                          # binning mechanism below is copied from scanpy's own
                          # _score_genes_bins, not just "inspired by" it (see comment there)

SEED = 42   # same value as 03_in_silico_perturb.py, recorded here explicitly for this script
random.seed(SEED)
np.random.seed(SEED)
rng = np.random.default_rng(SEED)
# ---------------------------------------------------------------------------


def load_perturb_module():
    spec = importlib.util.spec_from_file_location("perturb03", SCRIPT_03_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def get_tokenized_vocab(tok_dataset_path):
    """Every human Ensembl ID actually present in at least one cell's rank encoding in our
    tokenized dataset (not just anything in Geneformer's full ~25k-token vocabulary, most of
    which is never expressed in our data at all)."""
    from datasets import load_from_disk

    with open(TOKEN_DICT_PATH, "rb") as f:
        gene_token_dict = pickle.load(f)
    token_to_gene = {v: k for k, v in gene_token_dict.items() if k not in ("<pad>", "<mask>")}

    ds = load_from_disk(tok_dataset_path)
    seen_tokens = set()
    for row in ds["input_ids"]:
        seen_tokens.update(row)

    return {token_to_gene[t] for t in seen_tokens if t in token_to_gene}


def build_expression_matched_pool(candidates, tok_dataset_path, n_genes):
    """Bin genes by mean expression in our own FACS data, using the exact binning mechanism
    scanpy's sc.tl.score_genes uses internally (scanpy/tools/_score_genes.py,
    _score_genes_bins): rank genes by average expression, bin = rank // (n_pool / (n_bins-1)).
    Then draw n_genes total, split across bins in proportion to how many of the candidate
    genes fall in each bin, so the null sample's expression profile matches the candidates'.
    """
    adata = sc.read_h5ad(PREPPED_H5AD)
    mean_expr = np.asarray(adata.X.mean(axis=0)).ravel()
    expr_df = (
        pd.DataFrame({"ensembl_id": adata.var["ensembl_id"].values, "mean_expr": mean_expr})
        .dropna(subset=["ensembl_id"])
        .drop_duplicates(subset=["ensembl_id"])
    )

    vocab = get_tokenized_vocab(tok_dataset_path)
    samplable_pool = set(expr_df.loc[expr_df["ensembl_id"].isin(vocab), "ensembl_id"])
    expr_df = expr_df[expr_df["ensembl_id"].isin(vocab)]
    print(f"Gene pool for expression matching: {len(expr_df)} genes "
          f"(ortholog-mapped AND actually present in the tokenized dataset)")

    # add back any candidate genes missing from the samplable pool (e.g. OCT4/NANOG, see
    # MANUAL_MOUSE_SYMBOL above), using their true expression, purely so bin assignment for
    # those candidates is correct. Never added to samplable_pool, so never draw-eligible.
    missing_candidates = set(candidates) - set(expr_df["ensembl_id"])
    if missing_candidates:
        raw_adata = sc.read_h5ad(RAW_FACS_H5AD)
        extra_rows = []
        for gene in missing_candidates:
            mouse_sym = MANUAL_MOUSE_SYMBOL.get(gene)
            if mouse_sym is not None and mouse_sym in raw_adata.var_names:
                col = raw_adata[:, mouse_sym].X
                col = col.toarray().ravel() if hasattr(col, "toarray") else np.asarray(col).ravel()
                extra_rows.append({"ensembl_id": gene, "mean_expr": col.mean()})
                print(f"  {gene} ({mouse_sym}) not in the samplable pool (one2many ortholog, "
                      f"PROGRESS.md 4.4a), true raw expression looked up directly: "
                      f"mean={col.mean():.6f} across {raw_adata.n_obs} cells. Placed in a bin "
                      f"for reference only, not draw-eligible.")
            else:
                print(f"WARNING: candidate {gene} has no known mouse-symbol lookup, cannot "
                      f"place it in an expression bin at all.")
        if extra_rows:
            expr_df = pd.concat([expr_df, pd.DataFrame(extra_rows)], ignore_index=True)

    obs_avg = expr_df.set_index("ensembl_id")["mean_expr"]
    n_items = int(round(len(obs_avg) / (N_EXPRESSION_BINS - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items  # same formula as scanpy's _score_genes_bins

    candidates_present = [g for g in candidates if g in obs_cut.index]
    missing = set(candidates) - set(candidates_present)
    if missing:
        print(f"WARNING: candidate gene(s) could not be placed in any expression bin at all, "
              f"excluded from bin-matching: {missing}")

    candidate_cuts = obs_cut.loc[candidates_present]
    print("Candidate gene expression bins (bin index, higher = more highly expressed):")
    print(candidate_cuts.sort_values().to_string())

    cut_counts = candidate_cuts.value_counts().sort_index()
    is_candidate = obs_cut.index.isin(candidates)
    # only draw from genes that are both unambiguous 1:1 orthologs AND actually tokenized
    # (samplable_pool), never the manually-placed reference-only rows added above
    is_samplable = obs_cut.index.isin(samplable_pool)

    sampled = []
    remaining = n_genes
    cuts_list = list(cut_counts.items())
    for i, (cut, cnt) in enumerate(cuts_list):
        is_last = i == len(cuts_list) - 1
        n_this_bin = remaining if is_last else round((cnt / cut_counts.sum()) * n_genes)
        n_this_bin = min(n_this_bin, remaining)

        bin_pool = obs_cut[(obs_cut == cut) & ~is_candidate & is_samplable].index.tolist()
        if len(bin_pool) < n_this_bin:
            print(f"WARNING: bin {cut} only has {len(bin_pool)} non-candidate genes, "
                  f"requested {n_this_bin}; sampling with replacement for this bin")
            chosen = list(rng.choice(bin_pool, size=n_this_bin, replace=True))
        else:
            chosen = list(rng.choice(bin_pool, size=n_this_bin, replace=False))
        sampled.extend(chosen)
        remaining -= n_this_bin

    return sampled


def run_null(n_genes):
    p03 = load_perturb_module()
    genes = build_expression_matched_pool(p03.TFS_TO_PERTURB, p03.TOK_DATASET, n_genes)
    print(f"\nSampled {len(genes)} expression-matched random genes (seed={SEED}):")
    print(genes)

    null_dir = Path(OUT_DIR) / f"{OUT_PREFIX}_n{n_genes}"
    state_embs_dict = p03.get_state_embs(null_dir)

    rows = []
    for i, gene in enumerate(genes):
        gene_dir = null_dir / gene
        t0 = time.time()
        stats_csv = p03.run_one_perturbation([gene], state_embs_dict, gene_dir, gene)
        elapsed = time.time() - t0
        df = pd.read_csv(stats_csv)
        shift = df["Shift_to_goal_end"].iloc[0]
        print(f"  [{i + 1}/{len(genes)}] {gene}: {elapsed:.1f}s, Shift_to_goal_end={shift}")
        rows.append({"human_ensembl_id": gene, "goal_state_shift": shift, "runtime_s": round(elapsed, 1)})

    null_df = pd.DataFrame(rows)
    out_path = null_dir / f"{OUT_PREFIX}_n{n_genes}_shifts.csv"
    null_df.to_csv(out_path, index=False)
    return null_df, out_path


if __name__ == "__main__":
    print(f"SEED={SEED}, N_RANDOM_GENES={N_RANDOM_GENES}, N_EXPRESSION_BINS={N_EXPRESSION_BINS}")

    t0 = time.time()
    null_df, out_path = run_null(N_RANDOM_GENES)
    elapsed = time.time() - t0

    print(f"\nnull distribution written to {out_path}")
    print(f"total runtime: {elapsed:.1f}s")
    shifts = null_df["goal_state_shift"]
    print(
        f"distribution of shifts: mean={shifts.mean():.10f}, sd={shifts.std(ddof=1):.10f}, "
        f"min={shifts.min():.10f}, max={shifts.max():.10f}"
    )
