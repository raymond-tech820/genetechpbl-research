# FACS Limb Muscle MSC Youth Score v3.1

Minimal GitHub package for two frozen TMS FACS Limb Muscle MSC Youth Score candidates and their completed training and external-comparison evidence.

## Folder Layout

- `M1_raw_all_gene_denominator/`: standalone primary model package.
- `M4_raw_all_gene_denominator_pi_0_90/`: standalone exploratory comparator package.
- `training_experiment/`: combined training, controlled-ablation, freeze, and formal 999-permutation report with figures and key tables.
- `cross_assay_comparison/`: three reports covering GSE176206, TMS Droplet, and final M1-M4 synthesis, with figures and key tables.

## Model Roles

M1 is the primary deployable reference. M4 is a post-ablation high-stringency sensitivity comparator. External evidence did not justify upgrading M4, so these roles remain frozen.

## Quick Validation

```bash
Rscript M1_raw_all_gene_denominator/validate_model_package.R
Rscript M4_raw_all_gene_denominator_pi_0_90/validate_model_package.R
```

Both scorers require genes-by-samples raw pseudobulk counts and use the input matrix column sums as the raw all-gene CPM denominator. Read the selected model directory's README before application.

## Training Data

Raw training data are intentionally excluded from this GitHub package. The model was trained on the diaphragm-excluded TMS FACS Limb Muscle MSC dataset: 815 cells, 22,966 genes, and 14 mice.

Training data URL: `TODO_ADD_GOOGLE_DRIVE_OR_DATA_PORTAL_URL`

## Evidence Summary

- Formal 999-permutation primary statistic: M1 empirical p = 0.023; exploratory M4 p = 0.010.
- GSE176206: both models passed coverage, preserved baseline Young > Aged direction in both independent controls, and showed negative median SOKM-control contrasts against both controls. Comparisons were unpaired because animal labels are nested within treatment arms.
- TMS Droplet: both models had complete coverage, AUC 1.0, negative all-age and old-only correlations, and high donor agreement. Paired-bootstrap differences between M4 and M1 included zero.

## Interpretation Boundary

The models quantify a frozen FACS-relative young-like MSC transcriptional state. The evidence supports non-random age-label association, GSE directional and perturbation consistency, and within-TMS cross-assay transportability. It does not establish a universal aging clock, technical independence, rejuvenation, treatment safety, or formal Identity/Risk coupling.

## Dependencies

Model application and package self-tests use base R only. Training dependencies are documented in the experiment report but are not needed for frozen scoring.
