# PROGRESS.md

Running log for the Geneformer perturbation work. This file is carried back and forth
between Claude Code (technical work) and chat (research, interpretation, team comms).

**How to use this file**

- In Claude Code: update sections 3, 4, 5, 6, 10 as work happens. Keep entries short and factual.
- In chat: paste this file to get help with interpretation, framing, next steps, and team messages.
  Chat updates sections 2, 7, 8, 9 and hands the file back.
- Do not delete history. Append. Old entries are the audit trail.

Last updated: (fill in)
Updated by: (Claude Code | chat)

---

## 1. One-line status

> Current state: repo scaffolded, data downloaded but not yet verified. Nothing blocking.

---

## 2. Pre-registered analysis plan (fill in BEFORE running anything)

Written down in advance so results cannot be metric-shopped after the fact.

| Item | Decision | Locked? |
|---|---|---|
| Primary metric | (e.g. rank of Yamanaka factors among all tested perturbations by goal-state shift) | no |
| Secondary metric | | no |
| Null | random gene perturbations, N = ???, expression-matched? | no |
| Baseline | | no |
| What counts as a positive result | | no |
| What counts as a negative result | | no |
| Aggregation unit | donor level, not cell level | yes |
| Seed | | no |

Note: a negative result is a publishable outcome for this project. The supervisor explicitly
praises honest negatives. Define in advance what would falsify the approach.

---

## 3. Environment and setup log

| Date | Step | Status | Notes |
|---|---|---|---|
| | conda env `genetech-gf` created | | |
| | Geneformer cloned and installed | | version / commit hash: |
| | `torch.cuda.is_available()` | | GPU: |
| | Model checkpoint chosen | | V1-10M / V2-104M / V2-316M |
| | Matching token dictionary + gene median file | | must match checkpoint version |
| | git repo initialised, `.gitignore` in place | | pushing manually, no automation |
| | First push (scaffold only, no results needed) | | repo URL: |

### .gitignore contents (large files must never be staged)

```text
data/raw/
data/prepped/
data/tokenized/
Geneformer/
*.h5ad
*.mtx
__pycache__/
.ipynb_checkpoints/
```

---

## 4. Data log

### 4.1 Files received from Zihan (Google Drive)

Stored in `data/raw/`, gitignored.

| File | Downloaded | Notes |
|---|---|---|
| `tms_droplet_limb_msc_raw_counts.h5ad` | | the main input |
| `README_geneformer_input.md` | | written for this pipeline specifically, READ FIRST |
| `tms_droplet_limb_msc_h5ad_validation.txt` | | may already answer the raw-counts question |
| `raw_data_export_summary.txt` | | |
| `tms_droplet_limb_msc_cell_metadata.csv` | | cross-check against `.obs` |
| `tms_droplet_limb_msc_gene_metadata.csv` | | cross-check against `.var` |
| `.mtx` / `cells.txt` / `genes.txt` | skipped | redundant fallback format, only needed if the h5ad is malformed |

### 4.2 Reference artifacts from teammates

Stored in `data/external/`, small enough to commit.

| File | Source | Purpose | Received |
|---|---|---|---|
| `limb_msc_general_youth_score_v1_signature.csv` | Zihan repo | 100 genes, weights, young-high / old-high modules | |
| `limb_msc_general_youth_score_v1_calibration.json` | Zihan repo | mu, sigma, M_young, M_old (needed to actually compute his score) | |
| `score_limb_msc_youth.R` | Zihan repo | authoritative definition of the scoring formula | |
| Risk gene lists per axis | Kei | 3 axes: dedifferentiation, inflammation/SASP, proliferation | promised Tuesday, using candidate lists meanwhile |

### 4.3 Verification checklist (do this before any processing)

| Check | Result | Notes |
|---|---|---|
| `.X` contains integers, not decimals | | decimals means already normalized, WRONG input for Geneformer |
| `.raw` exists? If so, which layer holds raw counts | | some pipelines stash log-normalized data in `.raw` |
| Gene ID namespace in `.var` | | mouse symbols or mouse Ensembl? determines ortholog mapping key |
| `.obs` columns present | | need: mouse id, age, sex, cell type, (condition if applicable) |
| Cell count / gene count | | expect roughly 9,649 cells |
| Cell type is limb muscle MSC only | | |
| Young / old cell counts | | expect imbalance, roughly 1,468 young vs 8,181 old |

### 4.4 Processing log

| Date | Step | Status | Notes |
|---|---|---|---|
| | Ortholog mapping run | | genes retained: ??? / ??? (??%) |
| | Tokenized | | output path, n cells |

### 4.5 Open data questions

- [ ] Is the team standardizing on FACS or Droplet? Currently holding Droplet.
- [ ] Are Kei's final risk gene lists in?

---

## 5. Runs log

One row per actual run. Include seeds and settings so runs are reproducible.

| Run ID | Date | Script | Model | Genes perturbed | max_ncells | Seed | Outcome | Output path |
|---|---|---|---|---|---|---|---|---|
| | | | | | | | | |

### Failures and fixes

Keep these. They are the most useful thing to hand back to chat, and they become the
troubleshooting section of the final writeup.

| Date | What broke | Root cause | Fix |
|---|---|---|---|
| | | | |

---

## 6. Results so far

### 6.1 Headline numbers

(fill in as they arrive, with the metric name and the aggregation unit stated explicitly)

### 6.2 Ranking table

| Factor | Human Ensembl ID | Predicted shift toward young | Rank | Percentile versus null |
|---|---|---|---|---|
| | | | | |

### 6.3 Null distribution

- N random perturbations:
- Where the Yamanaka factors sit relative to the null:

### 6.4 Baseline comparison

- Baseline method:
- Does Geneformer beat it:

---

## 7. Interpretation and open scientific questions

(chat fills this in, Claude Code reads it for context)

- [ ] Does Geneformer's ranking agree with the real-data finding that SOKM did not rejuvenate
      aged MSCs (measured contrast -0.118)?
- [ ] Does Geneformer propose any factor NOT already tested experimentally that looks promising?
      This is the main thing the scoring models cannot do on their own.
- [ ] Are the top-ranked factors biologically plausible, or are they artifacts
      (e.g. highly expressed housekeeping genes that move any cell)?
- [ ] How sensitive is the ranking to the choice of goal state definition (old to young)?
- [ ] Do Geneformer's top hits overlap with the young-high genes in Zihan's 100-gene signature?
      Note: low overlap is not automatically bad. Chew's comparison doc argues that disjoint gene
      panels recovering the same axis is actually strong evidence. But we should know which case we are in.
      Requires mapping Zihan's mouse symbols through the same ortholog table.

---

## 8. Team coordination state

| Person | Owns | What I need from them | Status |
|---|---|---|---|
| Zihan | Youth Score, TMS **Droplet**, pseudobulk logCPM, R parser | raw counts (RECEIVED), signature + calibration files | |
| Kaile | Youth Score, TMS **FACS**, multi-model (gene signature / elastic net / transformer), `.h5ad` raw-count Python scorer | nothing directly, but hers is the cell-level raw-count interface if ever needed | |
| Kei | Risk Score, 3 axes, MSigDB Hallmark + SenMayo | finalized gene list per axis | promised Tuesday |
| Mashiro | Integration, Youth-vs-Risk state space, Safe Zone plot | agreement on hand-off format for my ranking output | needs a call |
| Chew | Supervisor | has reviewed Zihan and Kaile, not me yet | |

**Correction logged:** the single-cell raw-count `.h5ad` scorer described in
`YOUTH_SCORE_MODEL_INPUT_OUTPUT.md` is **Kaile's FACS** work, not Zihan's.
Zihan's deliverable is pseudobulk logCPM in R on Droplet. I mixed these up in an earlier
group message and corrected it.

### Agreed interfaces

- My output to Mashiro: a **ranking table**, not an expression matrix. See CLAUDE.md section 5.
- Format to confirm with Mashiro: columns, gene ID namespace (human Ensembl versus mouse symbol),
  and whether he wants per-cell or aggregated values.

---

## 9. Next actions

Keep this list short. Three to five items maximum.

- [ ] Read `README_geneformer_input.md` and `h5ad_validation.txt` from Zihan before writing any code
- [ ] Run the verification checklist in 4.3 and record the answers
- [ ] Set up git with the `.gitignore` above, push the scaffold early (supervisor weights runnability heavily)
- [ ] Fill in the pre-registered plan in section 2 before the first real run
- [ ] Get V1-10M running end to end on a small `max_ncells` to prove the pipeline

---

## 10. Session notes

Append a short entry per working session. Date, what was attempted, what happened, what is next.

### (date) session 1
- Attempted:
- Result:
- Next:
