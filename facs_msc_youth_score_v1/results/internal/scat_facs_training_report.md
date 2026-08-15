# Youth Score Training Report: scat_facs

Generated: 2026-07-17T09:50:59.111183+00:00

## Dataset

- Role: `primary`
- Modality: `facs`
- Tissue: `SCAT`
- Cell type: `mesenchymal stem cell of adipose`
- Young ages: `3m`
- Old ages: `18m, 24m`
- Cells after matrix-derived QC: `1570`
- Donors: `15`
- Young cells/donors: `737` / `7`
- Old cells/donors: `833` / `8`

## Model Selection

- Official model: `elastic_net`
- Evidence status: `internally_supported`
- Rule: Highest donor OOF ROC-AUC; within 0.02 use lower Brier; remaining ties prefer simpler models.
- Technical-only AUC: `0.6964`
- Official-model AUC: `0.9464`

## Donor-level Metrics

| model            |   roc_auc |   pr_auc |   balanced_accuracy |     f1 |   brier |   roc_auc_ci_low |   roc_auc_ci_high |   age_spearman |   donors |
|:-----------------|----------:|---------:|--------------------:|-------:|--------:|-----------------:|------------------:|---------------:|---------:|
| elastic_net      |    0.9464 |   0.9571 |              0.8036 | 0.8    |  0.0883 |           0.8214 |            1      |        -0.8356 |       15 |
| gene_transformer |    0.8571 |   0.8677 |              0.6161 | 0.6667 |  0.1937 |           0.625  |            1      |        -0.6684 |       15 |
| gene_signature   |    0.7857 |   0.7932 |              0.8036 | 0.8    |  0.1995 |           0.5179 |            1      |        -0.6454 |       15 |
| technical_only   |    0.6964 |   0.7385 |              0.6696 | 0.6667 |  0.2399 |           0.3929 |            0.9464 |        -0.4168 |       15 |

## Interpretation

The score is a cross-validated estimate of similarity to the TMS 3-month reference state. A higher score means more young-like expression within this tissue/cell-type model. It is not a safety score, causal rejuvenation measurement, clinical age, or treatment recommendation.

The primary unit of validation is the mouse, not the cell. Sex, sequencing plate, library complexity, and age are partially confounded in TMS. The `confound_limited` status means the model should remain exploratory even if its discrimination metric is high.

## Files

- `oof_cell_scores.parquet`: out-of-fold cell scores.
- `oof_donor_scores.csv`: donor-aggregated scores.
- `model_metrics.csv`: baseline and Transformer comparison.
- `selection.json`: selected model and confounding status.
- `folds/`: fold-specific tokenizers, calibration, baselines, and Transformer checkpoints.
- `figures/`: report plots.
