# FACS Limb Muscle MSC Youth Score v2.1 Model Card

## Intended Use

This is a frozen TMS FACS Limb Muscle MSC cohort-state score trained after excluding old-only `Muscle Diaphragm` cells. It is intended for donor/sample-level scoring of raw-count Limb Muscle MSC data after pseudobulk aggregation.

Primary model: `factorial_medium_original`.

Do not describe this model as a validated aging clock, a technical-confounding-free model, or an externally validated universal MSC aging score.

## Training Data

- Source: TMS FACS Limb Muscle MSC raw-count subset.
- Diaphragm cells excluded: 120.
- Final training cells: 815.
- Genes: 22,966.
- Donors/mice: 14.
- Age groups: 3m young, 18m/24m old.
- Young cells: 264.
- Old cells: 551.
- Diaphragm cells remaining: 0.
- Remaining subtissue standardized as `Forelimb_Hindlimb`.

## Frozen Training Procedure

1. Aggregate single-cell raw counts to mouse-level pseudobulk.
2. Filter genes at CPM > 1 in at least 2 mice.
3. Use edgeR TMM and limma voom for DE/ranking.
4. Fit factorial cell-means design `~ 0 + sex:age_group` plus age-only DE.
5. Compute female, male, common, and interaction age contrasts.
6. Exclude listed sex-linked genes.
7. Require female/male direction concordance and age-only/common direction concordance.
8. Rank candidates by `abs(common_logFC) * abs(common_t) * interaction_penalty`.
9. Select top 50 young-high and top 50 old-high genes for the primary Medium model.
10. Freeze full-data log2(CPM + 1) training mean/sd for each signature gene.
11. Calibrate raw module difference with training medians.

Training-set medians are fixed by calibration; apparent full-data separation is not validation evidence.

## Scoring Formula

```text
E_g = log2(CPM_g + 1)
Z_g = (E_g - mu_g) / s_g
M_young = weighted_mean(Z_g for young-high genes)
M_old   = weighted_mean(Z_g for old-high genes)
S_raw   = M_young - M_old
S_cal   = (S_raw - C_old) / (C_young - C_old)
```

Higher calibrated scores indicate a more young-like FACS training-cohort state.

## Key Internal Validation

Primary nested LOMO model, `factorial_medium_original`:

- AUC young vs old: 0.979.
- all-age Spearman: -0.633.
- old-only Spearman: 0.655.
- library-size Spearman: -0.521.
- zero-signature folds: 0.

Donor bootstrap for the primary model:

- AUC 95% CI: [0.878, 1.000].
- all-age rho 95% CI: [-0.838, -0.120].
- library rho 95% CI: [-0.794, 0.100].
- young-old median difference 95% CI: [0.155, 0.395].

Formal 999 nested permutation for the primary model:

- primary `abs_all_age_rho`: observed 0.633, empirical p = 0.153.
- supporting `abs(AUC - 0.5)`: observed 0.479, empirical p = 0.039.
- recorded young-minus-old median: empirical p = 0.002.
- valid null permutations: 999/999.
- failed permutations: 0.
- rank-deficient folds: 0.
- zero-signature folds: 0.

## Interpretation

The model provides reproducible evidence for binary young-old transcriptional state separation, but the preregistered primary permutation test does not support a robust continuous all-age trajectory beyond the practical randomized-label null.

Interpret as: frozen FACS young-old cohort-state separator.

Do not interpret as: validated continuous aging trajectory model or aging clock.

## Limitations

- Technical independence has not been demonstrated.
- The cohort has small donor count and weak-support factorial folds.
- The model does not establish monotonic 3m > 18m > 24m trajectory behavior.
- External validity has not been established.
- Application to other assays or datasets should be reported as transportability/sensitivity unless independently validated.

## Parser Validation

Exported parser reproduced training implementation at machine precision:

- minimum gene coverage: 1.0.
- minimum weighted coverage: 1.0.
- maximum calibrated-score difference: approximately 1e-15.
