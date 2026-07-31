# Geneformer perturbation arm

Zero-shot in-silico perturbation prediction for the cellular rejuvenation
project, using Geneformer-V1-10M on the TMS FACS limb-muscle MSC cohort.

Contributor: Jia Qi Choy

## What this arm does

Predicts how overexpressing candidate transcription factors shifts aged
limb-muscle MSCs toward a young reference state, and validates whether that
reference axis is biologically meaningful before any prediction is trusted.

## Key finding

The Old-to-Young embedding axis is confounded with sequencing depth
(detected-gene count, Pearson r = 0.69), and depth is inseparable from age in
this cohort, replicated across the FACS and Droplet assays. Fixed-length input
truncation reduces the confound but floors near r = 0.5 and does not survive
donor-level aggregation. The perturbation ranking is therefore not
interpretable on this cohort. The arm's contribution is this characterised
negative and a validation criterion: a usable aging axis must show donor-level
separation that survives conditioning on detected-gene count.

Full methods and results: [`work_progress.md`](./work_progress.md).

## Layout

- `src/` — pipeline (ortholog mapping, tokenization, perturbation) and analysis
  scripts (axis validation, deconfounding, token-length sweep)
- `work_progress.md` — full methods log and results
- `results/` — summary tables and figures (large artifacts gitignored)
- `data/external/` — teammate reference files
- `environment.yml` — conda environment

## Data

- Aging reference: Tabula Muris Senis (TMS) FACS limb-muscle MSC
- Reprogramming data: GSE176206 (Roux et al. 2022, *Cell Systems*)
- Model: Geneformer-V1-10M, human-pretrained (mouse genes ortholog-mapped)

---

# Setup notes

## Model choice

Geneformer ships several checkpoints. This arm uses the smallest, V1-10M, which
runs on a laptop GPU; larger checkpoints need a server.

| Model | Params | Input size | Rough VRAM | Use when |
|---|---|---|---|---|
| Geneformer-V1-10M | 10M | 2048 | ~2-4 GB | laptop / first pass (used here) |
| Geneformer-V2-104M | 104M | 4096 | ~10-16 GB | server GPU |
| Geneformer-V2-316M | 316M | 4096 | 24 GB+ | large server GPU only |

V1 and V2 use **different token dictionaries and gene-median files**. If the
model is switched, the matching dictionary files must be switched too; the
tokenizer takes them as arguments. Do not mix a V1 dictionary with a V2 model.
This arm uses V1-10M with the gc30M dictionary throughout.

## Environment

```bash
conda env create -f environment.yml
conda activate genetech-gf

# Geneformer installs from Hugging Face, not PyPI (needs git-lfs)
git lfs install
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer && pip install . && cd ..
```

Sanity check:

```bash
python -c "import geneformer, scanpy, torch; print('ok', torch.cuda.is_available())"
```

If `torch.cuda.is_available()` prints `False`, fix CUDA/PyTorch before running;
Geneformer needs a GPU to be usable. Note: transformers must be pinned below
version 5 (the environment file pins this); newer versions remove a symbol
Geneformer imports.

## Run order

```bash
python src/01_map_orthologs.py       # build mouse -> human ortholog table (once)
python src/02_prepare_tokenize.py    # map to human genes, tokenize (model_version="V1")
python src/03_in_silico_perturb.py   # in-silico overexpression + stats
```

Analysis scripts (`pseudobulk_axis.py`, `fixed_length_axis.py`,
`fixed_length_axis_cellset_control.py`, `token_length_sweep.py`,
`plot_figures_confound.py`) are documented step by step in `work_progress.md`.

## Notes

1. **Species.** Geneformer is human-pretrained; the data is mouse. Genes are
   mapped to human Ensembl IDs before tokenizing (step 1). One-to-one orthologs
   only are kept: 62.3% of genes, 68.9% of count mass. This species mismatch is
   central to the arm's finding, since the rank encoding uses human gene medians.

2. **Perturbation targets.** Yamanaka and related factors (Pou5f1/OCT4, Sox2,
   Klf4, Myc, plus Nanog, Lin28a, Myod1) are converted to human Ensembl IDs and
   perturbed with `overexpress`, not delete.

3. **Goal state.** The perturber measures shift toward a target cell state:
   start = old cells, goal = young cells, from the age labels. Positive
   `Shift_to_goal_end` means movement toward young. See `cell_states_to_model`
   in step 3.

4. **Dictionary/version match.** The tokenizer and perturber must be run with
   `model_version="V1"`; the default is V2, which silently mismatches the V1
   checkpoint and mis-maps gene tokens without erroring.

5. **Zero-shot.** The model is run zero-shot, not fine-tuned, which avoids
   fitting parameters to a small donor cohort.
