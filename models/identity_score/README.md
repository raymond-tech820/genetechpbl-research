# Identity Score v1 - code release

This directory contains the upload-ready implementation of the frozen mouse
MSC Identity Score v1 used for the 2026-08-06 GSE176206 handoff.

## What the score means

Identity Score v1 is a knowledge-driven expression-module score, not a trained
classifier. Higher values indicate stronger expression of the frozen mouse MSC
identity program relative to expression-matched background genes. It is not a
measure of clinical safety, MSC potency, causal rejuvenation, or tumour risk.

The frozen positive genes are:

`Cd34, Cd44, Cxcl12, Eng, Itgav, Lepr, Ly6a, Nt5e, Pdgfra, Prrx1, Thy1`.

`Vcam1` was reviewed but excluded from the v1 core.

## Frozen cell-level algorithm

1. Read non-negative raw counts from `layers/counts` and mouse gene symbols
   from `var/gene_name`.
2. Normalize each cell to 10,000 total counts and apply `log1p`.
3. Estimate mean normalized expression using known-animal control cells.
4. Divide all input genes into 24 expression bins.
5. For every identity gene, sample one eligible gene from the same bin in each
   of 100 deterministic background sets (seed `20260729`).
6. Primary score = mean expression of observed identity genes minus the mean
   score of the 100 matched background sets.
7. Rank sensitivity score = mean within-cell expression percentile of observed
   identity genes; undetected genes contribute zero.
8. Require at least 70% identity-gene coverage. The frozen GSE176206 run had
   11/11 genes and no missing scores.

The complete all-module frozen table is included because the original scorer
excluded every frozen module gene from the eligible background pool. Uploading
only the 11 identity genes would not reproduce the primary score exactly.

## Correct statistical unit

Scores are calculated for cells and aggregated using the median within:

```text
age_group x exact_treatment_arm x animal_label
```

Animal labels `1`, `2`, and `3` are reused across GSE176206 treatment arms.
They are nested labels, not global subject IDs. Cross-arm treatment comparisons
must therefore be unpaired. The code deliberately preserves `Tg+/Dox-` and
`Tg-/Dox+` as separate control arms and contains no paired-treatment function.

## Installation and run

From the `msc_identity_score_v1` directory:

```text
python -m pip install -e .
identity-score-v1 run --config config/gse176206.toml --force
```

The example configuration expects the input at
`data/GSE176206_msc_sokm.h5ad`. The raw dataset is not distributed in this
repository. You may instead provide any local path with `--input-h5ad`.

Optional portable path overrides:

```text
identity-score-v1 run --config config/gse176206.toml \
  --input-h5ad /path/to/GSE176206_msc_sokm.h5ad \
  --output-directory /path/to/identity_results \
  --force
```

## Main outputs

- `identity_cell_scores_minimal.parquet`: cell-level primary and rank scores.
- `identity_donor_by_condition.csv`: formal animal-level table.
- `identity_condition_summary.csv`: descriptive exact-arm summaries.
- `identity_unknown_cell_summary.csv`: unknown-animal cells, descriptive only.
- `identity_gene_set_qc.csv`: coverage and scoring QC.
- `validation_summary.json`: run-level assertions.

For formal integration, use `identity_primary_median` as the primary animal
summary, `identity_rank_median` as the sensitivity summary, and
`biological_unit_id` as the join key.

## Tests

```powershell
python -m pytest
```

The tests explicitly verify that reused animal labels in different treatment
arms remain distinct, unpaired biological units.


To compare a completed run with the frozen small reference tables:

```text
python scripts/verify_release.py --candidate-directory /path/to/identity_results
```

For a full cell-level regression check, also supply the 2026-08-06 minimal
cell-score parquet with `--reference-cell-scores`.