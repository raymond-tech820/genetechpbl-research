# Geneformer: MSC in-silico perturbation

My piece of the cellular-rejuvenation project: use **Geneformer** to predict how
transcription factors / reprogramming factors shift limb-muscle MSCs toward a
"younger" state, and flag when they also push toward the risky (identity-loss /
proliferative) state.

- **Cell type:** limb muscle–derived MSCs (same population as the Youth Score work)
- **Aging reference:** Tabula Muris Senis (TMS)
- **Reprogramming data:** GSE176206 (Gill et al. 2022, *Cell Systems*, Dox pulse/chase)
- **Model:** Geneformer (human-pretrained → mouse genes need ortholog mapping)

---

## 0. Pick the right model for your GPU

Geneformer ships several checkpoints. The default is large; on a small GPU use the
smallest one first and only scale up on a server.

| Model | Params | Input size | Rough VRAM | Use when |
|---|---|---|---|---|
| Geneformer-V1-10M | 10M | 2048 | ~2–4 GB | laptop / first pass (start here) |
| Geneformer-V2-104M | 104M | 4096 | ~10–16 GB | server GPU |
| Geneformer-V2-316M | 316M (default) | 4096 | 24 GB+ | big server GPU only |

Start with **V1-10M** to get the whole pipeline working end-to-end, then swap the
model path to a bigger one on a server if you have time. Getting a result you can
explain beats getting the biggest model to load.

> Note: V1 and V2 use **different token dictionaries and gene-median files**. If you
> switch models, you must switch the matching dictionary files too (the tokenizer
> takes them as arguments). Don't mix a V1 dictionary with a V2 model.

---

## 1. Environment

```bash
# clone this project, then:
conda env create -f environment.yml
conda activate geneformer-msc

# install Geneformer itself (it lives on Hugging Face, not PyPI)
# needs git-lfs: https://git-lfs.com
git lfs install
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
pip install .
cd ..
```

Sanity check:

```bash
python -c "import geneformer, scanpy, torch; print('ok', torch.cuda.is_available())"
```

If `torch.cuda.is_available()` prints `False`, fix CUDA/PyTorch before going further,
Geneformer needs a GPU to be usable.

---

## 2. Run order

```bash
# 1. build the mouse -> human ortholog table (run once, caches to a CSV)
python src/01_map_orthologs.py

# 2. load data, attach ensembl_id + n_counts, map to human genes, tokenize
python src/02_prepare_tokenize.py

# 3. run in-silico perturbation of candidate TFs and get the stats
python src/03_in_silico_perturb.py
```

Each script has a `# TODO` block at the top for the paths / gene lists you need to fill in.

---

## 3. The gotchas (read before you start)

1. **Species.** Geneformer is human-pretrained; TMS and GSE176206 are mouse. Everything
   must be converted to **human Ensembl gene IDs** before tokenizing. That's what
   step 1 is for. Genes with no 1:1 human ortholog get dropped, expect to lose some.

2. **Yamanaka factors are the perturbation targets.** Convert *Pou5f1/Oct4, Sox2, Klf4,
   Myc* (+ any TFs the biology lead flags) to their **human Ensembl IDs** and feed those
   as the genes to perturb. Perturb `overexpress` (reprogramming = forcing them on), not delete.

3. **"Younger" as a goal state.** Geneformer's perturber can measure the shift *toward a
   target cell state*. Define start = old cells, goal = young cells (from the age labels),
   and it scores each TF by how much it moves cells toward "young." That's our Youth axis
   directly. (See the `cell_states_to_model` argument in step 3.)

4. **API drift.** Geneformer's function arguments have changed across versions. The three
   scripts follow the documented pattern, but **before running, open the example notebook
   in the version you cloned** (`Geneformer/examples/in_silico_perturbation.ipynb` and
   `.../tokenizing_scRNAseq_data.ipynb`) and match the argument names. Treat my scripts as
   a skeleton, the notebook as ground truth.

5. **Zero-shot first.** You do **not** need to fine-tune to get first results, run the
   pretrained model zero-shot. Only fine-tune a cell-state classifier later if the
   zero-shot shift isn't clean enough.

---

## 4. What "done" looks like for this piece

A table: for each candidate TF, its predicted shift toward "young" (Youth axis) and toward
the risky state (Risk axis). The winners are the ones that move cells up on Youth without
moving them up on Risk, i.e. the perturbations that land in the Safe Zone.
