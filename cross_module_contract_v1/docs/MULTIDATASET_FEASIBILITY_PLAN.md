# Multi-Dataset Feasibility Plan

No new data are downloaded or merged in this phase. This framework classifies future MSC datasets as `eligible_training`, `eligible_internal_validation`, `eligible_external_validation`, `eligible_cross_assay_sensitivity`, `descriptive_only`, `not_eligible`, or `TBD`.

## Intake checks

Record independent donor count, donors per age group, sex balance, age-sex confounding, dataset-age confounding, assay-age confounding, tissue/subtissue compatibility, MSC annotation compatibility, raw-count availability, donor-resolved and paired-sample metadata, genome build, gene identifier namespace, transcript complexity by age, Kaile feature coverage, Zihan signature coverage, and Kei Risk coverage.

## Analysis rules

- No simple concatenation.
- Use donor-level inference and within-dataset age contrasts.
- Preserve dataset labels and use leave-one-dataset-out validation.
- Use dataset-aware covariates and an explicit transcript-complexity audit.
- Report dataset-specific effects; more cells do not equal more independent donors.
- Harmony-corrected values cannot be count-model input.
- Do not claim multi-dataset harmonization is complete until these checks and validation are run.

