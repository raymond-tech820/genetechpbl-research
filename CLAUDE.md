# CLAUDE.md

Project context for Claude Code. Read this first in every session.

## Style rules (apply to all output, including code comments and markdown)

- Never use em dashes. Use commas, colons, parentheses, or separate sentences.
- Plain language in documentation. The team includes people with no ML background and people with no biology background.
- Be honest about negative results. This project's supervisor explicitly praises honest negatives and penalizes overselling.

---

## 1. What this project is

A 6-week team PBL research project: **AI for Safe Cellular Rejuvenation**.

Core question: can we identify reprogramming conditions that make aged cells
transcriptionally younger **without** activating cancer-like or unsafe cell states?

The team's framework has two axes and one predictor:

- **Youth Score** (axis 1): how young does this cell look? Built by teammates.
- **Risk Score** (axis 2): how much does it look dedifferentiated / inflamed / hyperproliferative? Built by a teammate.
- **Geneformer perturbation prediction**: THIS REPO. Which transcription factors, if
  overexpressed, are predicted to move cells toward "young"?

The target is the "Safe Zone": interventions that raise Youth without raising Risk.

## 2. My specific role (owner of this repo)

I own the **prediction** side, not the scoring side. Nobody else on the team is doing this.

- Input: mouse limb-muscle MSC single-cell data (raw counts) plus a list of candidate genes.
- Method: **Geneformer, pretrained, run zero-shot**. No fine-tuning, no training on our mice.
- Output: a **ranking table** of candidate factors by predicted shift toward the young cell state,
  with a null distribution and a simple baseline for comparison.

### Why zero-shot specifically

A teammate trained a gene-token transformer from scratch on this data and it underperformed
simple models (AUC 0.857 / 0.875 versus 1.000 for a gene signature). The supervisor's review
noted a transformer is oversized for ~15 donors. Zero-shot with a model pretrained on ~100M
cells avoids that failure mode entirely, because we never train on our 12 to 15 mice.
This is the main scientific justification for the approach. Keep it.

## 3. Data

| Dataset | What it is | Role here |
|---|---|---|
| TMS (Tabula Muris Senis) Droplet limb-muscle MSC | Mouse aging atlas, 12 mice, ~9,649 cells, raw counts, shared by teammate Zihan | Defines the old versus young reference states |
| GSE176206 (Gill et al. 2022, Cell Systems) | Partial reprogramming pulse/chase, Young/Aged x Control/SOKM | The reprogramming arm, used for comparison against real measured effects |

Age definition used by the team: Young = 3 months. Old = 18 / 21 / 24 months.

### Known data caveats (do not paper over these)

- TMS Droplet limb MSC has **only 2 young mice and both are female**, versus 10 old mice of
  mixed sex. Any age effect is estimated within young females. Do not claim sex independence.
- Young/Old cell counts are very imbalanced (1,468 young versus 8,181 old).
- There is an unresolved team decision about whether **FACS** or **Droplet** is the primary
  assay. The raw data I currently have is Droplet. Check before assuming.

## 4. Critical technical constraints

### 4.1 Geneformer needs RAW counts

Geneformer builds its own rank-value encoding internally. Never feed it normalized data.
If a matrix contains decimals it is already normalized and is the wrong input.
Raw counts are integers. Also check `adata.raw.X` versus `adata.X`: some pipelines store
log-normalized values in `.raw`, which would silently double-normalize.

### 4.2 Species mismatch: Geneformer is human-pretrained, our data is mouse

All genes must be mapped to **human Ensembl IDs** before tokenizing. Use one-to-one orthologs
only, at least to start. Expect to lose some genes. Report the retention rate.

### 4.3 Model size versus GPU

Default Geneformer checkpoint is V2-316M and will not fit a small GPU.
Start with **V1-10M** (input size 2048) to get the pipeline working end to end,
then scale up on a server if time allows.
V1 and V2 use **different token dictionaries and gene-median files**. Never mix a V1
dictionary with a V2 model.

### 4.4 Geneformer API drift

Argument names for `TranscriptomeTokenizer`, `InSilicoPerturber`, and
`InSilicoPerturberStats` have changed between versions. Always cross-check against the
example notebooks in the cloned version:
`Geneformer/examples/tokenizing_scRNAseq_data.ipynb` and
`Geneformer/examples/in_silico_perturbation.ipynb`.
Treat the notebooks as ground truth and this repo's scripts as a skeleton.

### 4.5 Perturbation direction

Reprogramming means forcing factors ON. Use `perturb_type="overexpress"`, not delete.
Yamanaka factors: Oct4 (POU5F1), Sox2, Klf4, Myc. Convert to human Ensembl IDs.

## 5. The interface seam with the rest of the team

Geneformer's in-silico perturbation outputs an **embedding-space shift toward a goal state**.
It does NOT output a reconstructed expression matrix. Therefore its output **cannot be fed
directly into the teammates' Youth Score or Risk Score parsers**, which expect expression data.

Do not attempt to reconstruct an expression matrix from the embedding shift. That path is
fragile and was explicitly rejected. The agreed design instead:

1. A teammate (Mashiro) scores the **real** GSE176206 cells with Youth + Risk. No Geneformer involved.
2. This repo produces an **independent ranking** of factors by predicted shift toward young.
3. We compare the two rankings.

Useful cross-check: a teammate found SOKM did **not** move aged MSCs toward youth
(contrast -0.118). If Geneformer independently also predicts SOKM is a weak mover,
that is two methods agreeing. That is a real result worth reporting.

## 6. Standards the supervisor enforces

These came from written review feedback on teammates' work. Build them in from the start.

1. **Ship runnable code, not markdown.** The harshest criticism given to a teammate was that
   their deliverable had no code, no weights, no example input, no seed. Push the repo early.
2. **Baselines are mandatory.** A strong-looking number with no reference point is not
   interpretable. Always include a simple baseline.
3. **Nulls need enough iterations.** A teammate ran 20 permutations, so the minimum attainable
   p-value was 0.048, which is architectural rather than evidentiary. Use enough draws.
4. **Aggregate to donor level.** Cells are not independent biological replicates.
   Thousands of cells do not replace 12 mice.
5. **Pre-register the primary metric** in the README before running, to avoid metric-shopping.
6. **Set seeds everywhere** randomness enters.

## 7. Repo layout

```text
/home/user/MIT-PBL/genetech
  CLAUDE.md                  this file
  PROGRESS.md                running log, synced back and forth with chat
  README.md                  setup, run order, pre-registered metric
  environment.yml
  data/
    raw/                     teammate-provided .h5ad and metadata (gitignored, large)
    external/                Kei's risk gene lists, Zihan's youth signature
    mouse_to_human_orthologs.csv
    prepped/                 cleaned .h5ad ready for tokenizing
    tokenized/               Geneformer tokenizer output
  src/
    01_map_orthologs.py
    02_prepare_tokenize.py
    03_in_silico_perturb.py
    04_null_distribution.py  random-gene null
    05_baseline.py           simple non-transformer baseline
    06_rank_and_report.py    final ranking table for hand-off
  results/
  notebooks/
```

## 8. Environment

- WSL, conda env `genetech-gf`, working dir `/home/user/MIT-PBL/genetech`
- Geneformer installs from Hugging Face, not PyPI. Needs `git-lfs`.
  `git lfs install && git clone https://huggingface.co/ctheodoris/Geneformer && cd Geneformer && pip install .`
- I have an RTX 4060 laptop GPU. Plan around limited VRAM: use V1-10M, small
  `forward_batch_size`, and cap `max_ncells` for first runs.
- There is also a shared GPU server available for larger runs.

## 9. Working preferences

- I am a bachelor's student with a strong ML and systems background and **no biology background**.
  Explain biology in plain language. Do not assume I know gene names, cell types, or assay jargon.
- Explain what a step is for before writing the code for it.
- Flag when something is a project-specific judgement call versus a standard method,
  because I have to defend these choices to a supervisor and to biologist teammates.
- Prefer starting simple and escalating only if the simple version leaves something on the table.
