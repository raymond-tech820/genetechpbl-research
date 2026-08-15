# Youth Score Training Report: limb_facs

Generated: 2026-07-17T09:51:02.373638+00:00

## Dataset

- Role: `secondary`
- Modality: `facs`
- Tissue: `Limb_Muscle`
- Cell type: `mesenchymal stem cell`
- Young ages: `3m`
- Old ages: `18m, 24m`
- Cells after matrix-derived QC: `935`
- Donors: `14`
- Young cells/donors: `264` / `6`
- Old cells/donors: `671` / `8`

## Model Selection

- Official model: `gene_signature`
- Evidence status: `internally_supported`
- Rule: Highest donor OOF ROC-AUC; within 0.02 use lower Brier; remaining ties prefer simpler models.
- Technical-only AUC: `0.6875`
- Official-model AUC: `1.0000`

## Donor-level Metrics

| model            |   roc_auc |   pr_auc |   balanced_accuracy |     f1 |   brier |   roc_auc_ci_low |   roc_auc_ci_high |   age_spearman |   donors |
|:-----------------|----------:|---------:|--------------------:|-------:|--------:|-----------------:|------------------:|---------------:|---------:|
| gene_signature   |    1      |   1      |              1      | 1      |  0.0467 |           1      |                 1 |        -0.8626 |       14 |
| elastic_net      |    0.9792 |   0.9762 |              0.8542 | 0.8333 |  0.0584 |           0.8958 |                 1 |        -0.8815 |       14 |
| gene_transformer |    0.875  |   0.9167 |              0.9167 | 0.9091 |  0.1359 |           0.5833 |                 1 |        -0.6094 |       14 |
| technical_only   |    0.6875 |   0.6194 |              0.7917 | 0.7692 |  0.2427 |           0.3745 |                 1 |        -0.4735 |       14 |

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
