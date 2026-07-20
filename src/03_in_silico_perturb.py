"""
Step 3. In-silico perturbation of candidate TFs.

Overexpress each reprogramming factor / TF in silico and measure how far it
shifts cells toward the "young" cell state (our Youth axis). Repeat with a
"risky" goal state (proliferative / identity-loss) to get the Risk axis.

Run:     python src/03_in_silico_perturb.py
Output:  results/  (per-TF shift stats)

IMPORTANT: argument names for InSilicoPerturber / InSilicoPerturberStats and the
structure of `cell_states_to_model` have changed between Geneformer versions.
Open examples/in_silico_perturbation.ipynb in the version you cloned and match it.
This script is the skeleton, the notebook is ground truth.
"""

from pathlib import Path

# ---------------------------------------------------------------- TODO -----
MODEL_DIR = "Geneformer/gf-6L-30M-i2048"   # <- V1-10M start; swap for a V2 path on a server
TOK_DATASET = "data/tokenized/tms_droplet_limb_msc.dataset"
OUT_DIR = "results"
OUT_PREFIX = "tms_limb_msc_tf_perturb"

# Human Ensembl IDs of the factors to test (convert with the ortholog table).
# Yamanaka factors: Oct4/POU5F1, Sox2, Klf4, Myc  (overexpress = force them on)
TFS_TO_PERTURB = [
    # "ENSG00000204531",  # POU5F1 (OCT4)   <- fill in real IDs
    # "ENSG00000181449",  # SOX2
    # "ENSG00000136826",  # KLF4
    # "ENSG00000136997",  # MYC
]

# Define the states to model the shift toward. This uses an obs column we kept
# during tokenizing (e.g. "age") and the value that marks young vs old cells.
# Youth axis: start from old cells, measure shift toward the "young" state.
STATE_KEY = "age"
START_STATE = "old"     # value in the STATE_KEY column marking aged cells
GOAL_STATE = "young"    # value marking young cells
# ---------------------------------------------------------------------------


def run():
    from geneformer import InSilicoPerturber, InSilicoPerturberStats

    Path(OUT_DIR).mkdir(exist_ok=True)

    # cell_states_to_model: shift from START_STATE toward GOAL_STATE on STATE_KEY.
    # (Confirm exact dict shape in your version's example notebook.)
    cell_states_to_model = {
        "state_key": STATE_KEY,
        "start_state": START_STATE,
        "goal_state": GOAL_STATE,
        "alt_states": [],
    }

    isp = InSilicoPerturber(
        perturb_type="overexpress",       # forcing reprogramming factors ON
        genes_to_perturb=TFS_TO_PERTURB,
        model_type="Pretrained",          # zero-shot, no fine-tuning
        num_classes=0,
        emb_mode="cell",
        cell_states_to_model=cell_states_to_model,
        max_ncells=2000,                  # keep small for a first run / small GPU
        forward_batch_size=16,            # lower if you hit OOM on a 4GB GPU
        nproc=4,
    )
    isp.perturb_data(
        model_directory=MODEL_DIR,
        input_data_file=TOK_DATASET,
        output_directory=OUT_DIR,
        output_prefix=OUT_PREFIX,
    )

    # Summarize into a per-TF "shift toward goal state" table.
    stats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed=TFS_TO_PERTURB,
        cell_states_to_model=cell_states_to_model,
    )
    stats.get_stats(
        input_data_directory=OUT_DIR,
        null_dist_data_directory=None,
        output_directory=OUT_DIR,
        output_prefix=f"{OUT_PREFIX}_stats",
    )
    print(f"done, see {OUT_DIR}/{OUT_PREFIX}_stats*")


if __name__ == "__main__":
    if not TFS_TO_PERTURB:
        raise SystemExit("Fill in TFS_TO_PERTURB with human Ensembl IDs first (see TODO).")
    run()
