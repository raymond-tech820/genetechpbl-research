"""
Step 3. In-silico perturbation of candidate TFs.

Overexpress reprogramming factor(s) in silico and measure how far it shifts
cells toward the "young" cell state (our Youth axis). Repeat with a "risky"
goal state (proliferative / identity-loss) to get the Risk axis (not yet
implemented here).

Run:     python src/03_in_silico_perturb.py
Output:  results/  (per-TF or per-combo shift stats)

IMPORTANT: argument names for InSilicoPerturber / InSilicoPerturberStats / EmbExtractor
and the structure of `cell_states_to_model` have changed between Geneformer versions.
Ground truth is examples/in_silico_perturbation.ipynb in the cloned Geneformer version,
not this script. Findings from reading that notebook and the library source are recorded
in PROGRESS.md 4.4a/4.4b/5:
  - genes_to_perturb as a list is perturbed as ONE combined condition, not gene-by-gene.
    That's why this script has two modes (see MODE below).
  - InSilicoPerturber requires `state_embs_dict` whenever `cell_states_to_model` is set;
    it is NOT computed automatically. It's built here with EmbExtractor.get_state_embs().
  - model_version="V1" must be passed to EmbExtractor, InSilicoPerturber, AND
    InSilicoPerturberStats separately (each defaults to V2 / the gc104M dictionary on its
    own). Passing emb_mode="cell" without model_version="V1" is not enough by itself for
    InSilicoPerturber/EmbExtractor's *token dictionary* choice, only for emb_mode.
  - InSilicoPerturberStats.get_stats() scans its entire input_data_directory for every
    file matching pickle_suffix ("_raw.pickle" by default), it does not filter by
    output_prefix. Running more than one gene/combo into the same directory would silently
    mix their raw pickles together. Each run below therefore gets its own subdirectory.
"""

import random
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

# ---------------------------------------------------------------- TODO -----
MODEL_DIR = "Geneformer/Geneformer-V1-10M"   # actual checkpoint dir in the cloned repo (PROGRESS.md 3)
TOK_DATASET = "data/tokenized/tms_facs_limb_msc.dataset"
OUT_DIR = "results"
OUT_PREFIX = "tms_facs_limb_msc_tf_perturb"

# Reproducibility. NOTE: Geneformer's own cell-subsampling for max_ncells
# (perturber_utils.downsample_and_sort) hardcodes `shuffle(seed=42)` internally and is NOT
# controlled by this constant, it will use 42 regardless of what SEED is set to here. SEED
# covers everything else we control (numpy/torch/random calls in this script) and is set
# here mainly for the record, since the forward pass itself is deterministic anyway (the
# model is loaded in eval mode inside perturb_data, so its own dropout is disabled).
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

# "combo": one InSilicoPerturber call, all of TFS_TO_PERTURB overexpressed together as a
#          single combined condition (e.g. testing SOKM as a group, like the real GSE176206
#          experiment). Output: one goal-state-shift number for the whole combination.
# "individual": loop over TFS_TO_PERTURB, one InSilicoPerturber call per gene (each run in
#          its own output subdirectory, see module docstring), then assemble the per-gene
#          goal-state shifts into a single ranking table. This is the primary deliverable
#          format: a ranking table of candidate factors by predicted shift toward the young
#          cell state.
MODE = "individual"

# These two cell counts are deliberately separate constants, not one shared value, because
# they answer different questions:
#   - STATE_EMB_NCELLS controls how many real Old cells and real Young cells go into the
#     mean-embedding centroids that DEFINE "toward young" in embedding space (EmbExtractor).
#     Too few cells here makes the Old/Young centroids themselves noisy, which would corrupt
#     every downstream shift measurement regardless of which gene is perturbed. This should
#     use as many cells as feasible.
#   - PERTURB_NCELLS controls how many (already-Old) cells actually get the gene(s) forced on
#     and re-embedded, i.e. how many individual shift measurements get averaged into the
#     final number. This is mainly a compute-cost knob (each additional cell costs one more
#     forward pass), not a centroid-stability concern.
# Both are set to the full 815-cell tokenized dataset (551 Old, 264 Young after diaphragm
# exclusion, PROGRESS.md 4.3c/4.4) for this run; PERTURB_NCELLS is further capped by
# filter_data to the 551 Old cells regardless of the number set here.
STATE_EMB_NCELLS = 815
PERTURB_NCELLS = 815

# Human Ensembl IDs of the factors to test. See PROGRESS.md 4.4a for how these were chosen,
# including the OCT4/Nanog one2many ambiguity and how each was resolved or left unresolved.
TFS_TO_PERTURB = [
    "ENSG00000204531",  # POU5F1 (OCT4) -- resolved candidate, PROGRESS.md 4.4a (not ENSG00000212993/POU5F1B)
    "ENSG00000181449",  # SOX2
    "ENSG00000136826",  # KLF4
    "ENSG00000136997",  # MYC
    "ENSG00000111704",  # NANOG -- resolved candidate, PROGRESS.md 4.4a (not ENSG00000255192, absent from V1 vocab)
    "ENSG00000131914",  # LIN28A
    "ENSG00000129152",  # MYOD1
]

# Human-readable labels for the ranking table / reporting, keyed by the IDs above.
GENE_SYMBOLS = {
    "ENSG00000204531": "POU5F1",
    "ENSG00000181449": "SOX2",
    "ENSG00000136826": "KLF4",
    "ENSG00000136997": "MYC",
    "ENSG00000111704": "NANOG",
    "ENSG00000131914": "LIN28A",
    "ENSG00000129152": "MYOD1",
}

# Define the states to model the shift toward. This uses an obs column we kept during
# tokenizing (see KEEP_OBS in 02_prepare_tokenize.py) and the value that marks young vs old
# cells. Youth axis: start from old cells, measure shift toward the "young" state.
# Checked directly against the loaded FACS h5ad (PROGRESS.md 4.3b): the column with a
# young/old label is "age_group", with exact values "Old" / "Young" (capitalized).
# "age" is NOT this column, it holds "3m" / "18m" / "24m" instead.
STATE_KEY = "age_group"
START_STATE = "Old"     # value in the STATE_KEY column marking aged cells
GOAL_STATE = "Young"    # value marking young cells
CELL_STATES_TO_MODEL = {
    "state_key": STATE_KEY,
    "start_state": START_STATE,
    "goal_state": GOAL_STATE,
    "alt_states": [],
}
# ---------------------------------------------------------------------------


def get_state_embs(out_dir):
    """Mean embedding position of real Old cells and real Young cells, needed by
    InSilicoPerturber to know what "toward young" means in embedding space. Recomputed on
    every invocation for simplicity; if this script is used for many individual-gene runs
    later, caching this across runs would save time since it doesn't depend on which gene
    is perturbed."""
    from geneformer import EmbExtractor

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    embex = EmbExtractor(
        model_type="Pretrained",   # zero-shot: pretrained checkpoint only, never fine-tuned on our mice
        num_classes=0,
        emb_mode="cell",           # V1 has no <cls> token
        max_ncells=STATE_EMB_NCELLS,
        emb_layer=0,
        summary_stat="exact_mean",
        forward_batch_size=16,
        model_version="V1",       # loads the gc30M dictionary matching V1-10M (PROGRESS.md 4.4/4.4a)
        nproc=4,
    )
    return embex.get_state_embs(
        CELL_STATES_TO_MODEL,
        MODEL_DIR,
        TOK_DATASET,
        str(out_dir),
        "state_embs",
    )


def run_one_perturbation(genes, state_embs_dict, run_dir, run_prefix):
    """Run a single InSilicoPerturber + InSilicoPerturberStats pass for `genes` (perturbed
    together as one condition), scoped to its own `run_dir` so InSilicoPerturberStats's
    directory-wide pickle scan can't pick up another run's files. Returns the stats CSV path."""
    from geneformer import InSilicoPerturber, InSilicoPerturberStats

    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    isp = InSilicoPerturber(
        perturb_type="overexpress",       # forcing reprogramming factors ON
        genes_to_perturb=genes,
        model_type="Pretrained",          # zero-shot, no fine-tuning
        num_classes=0,
        emb_mode="cell",
        filter_data={STATE_KEY: [START_STATE]},   # only perturb cells actually starting in "Old"
        cell_states_to_model=CELL_STATES_TO_MODEL,
        state_embs_dict=state_embs_dict,
        max_ncells=PERTURB_NCELLS,
        emb_layer=0,
        forward_batch_size=16,             # lower if you hit OOM on a small GPU
        model_version="V1",                # gc30M dictionary, matches V1-10M
        nproc=4,
    )
    isp.perturb_data(
        model_directory=MODEL_DIR,
        input_data_file=TOK_DATASET,
        output_directory=str(run_dir),
        output_prefix=run_prefix,
    )

    stats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed=genes,
        cell_states_to_model=CELL_STATES_TO_MODEL,
        model_version="V1",                # gc30M dictionary, matches V1-10M (was missing before, PROGRESS.md 5)
    )
    stats.get_stats(
        input_data_directory=str(run_dir),
        null_dist_data_directory=None,
        output_directory=str(run_dir),
        output_prefix=f"{run_prefix}_stats",
    )
    return run_dir / f"{run_prefix}_stats.csv"


def run_combo():
    combo_dir = Path(OUT_DIR) / f"{OUT_PREFIX}_combo"
    state_embs_dict = get_state_embs(combo_dir)
    stats_csv = run_one_perturbation(TFS_TO_PERTURB, state_embs_dict, combo_dir, OUT_PREFIX)
    print(f"combo run done, see {stats_csv}")
    return stats_csv


def run_individual():
    individual_dir = Path(OUT_DIR) / f"{OUT_PREFIX}_individual"
    state_embs_dict = get_state_embs(individual_dir)

    rows = []
    for gene in TFS_TO_PERTURB:
        gene_dir = individual_dir / gene
        t_gene0 = time.time()
        stats_csv = run_one_perturbation([gene], state_embs_dict, gene_dir, gene)
        gene_elapsed = time.time() - t_gene0
        df = pd.read_csv(stats_csv)
        shift = df["Shift_to_goal_end"].iloc[0]
        print(f"  {GENE_SYMBOLS.get(gene, gene)} ({gene}): {gene_elapsed:.1f}s, Shift_to_goal_end={shift}")
        rows.append(
            {
                "gene_symbol": GENE_SYMBOLS.get(gene, gene),
                "human_ensembl_id": gene,
                "goal_state_shift": shift,
                "runtime_s": round(gene_elapsed, 1),
            }
        )

    ranking = (
        pd.DataFrame(rows)
        .sort_values("goal_state_shift", ascending=False)
        .reset_index(drop=True)
    )
    ranking.insert(0, "rank", ranking.index + 1)
    ranking = ranking[["rank", "gene_symbol", "human_ensembl_id", "goal_state_shift", "runtime_s"]]

    out_path = individual_dir / f"{OUT_PREFIX}_ranking.csv"
    ranking.to_csv(out_path, index=False)
    print(f"individual mode done, ranking table written to {out_path}")
    print(ranking.to_string(index=False))
    return out_path


if __name__ == "__main__":
    if not TFS_TO_PERTURB:
        raise SystemExit("Fill in TFS_TO_PERTURB with human Ensembl IDs first (see TODO).")

    Path(OUT_DIR).mkdir(exist_ok=True)
    print(f"SEED={SEED} (see comment above: this does not override Geneformer's own "
          f"hardcoded seed=42 for max_ncells subsampling, recorded for the record regardless)")

    has_cuda = torch.cuda.is_available()
    if has_cuda:
        torch.cuda.reset_peak_memory_stats()
    t0 = time.time()

    if MODE == "combo":
        run_combo()
    elif MODE == "individual":
        run_individual()
    else:
        raise SystemExit(f"Unknown MODE: {MODE!r}, expected 'combo' or 'individual'.")

    elapsed = time.time() - t0
    print(f"total runtime: {elapsed:.1f}s")
    if has_cuda:
        peak_gb = torch.cuda.max_memory_allocated() / (1024 ** 3)
        print(f"peak CUDA memory allocated: {peak_gb:.2f} GB")
    else:
        print("CUDA not available, ran on CPU (no VRAM figure to report)")
