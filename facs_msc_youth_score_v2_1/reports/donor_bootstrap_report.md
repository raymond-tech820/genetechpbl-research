# FACS Youth Score v2.1 Step 16: Donor Bootstrap Uncertainty

Bootstrap iterations: 10000 donor-level resamples with replacement per model.

Purpose: quantify uncertainty for predefined candidate/baseline models. This analysis is not a winner-selection or technical-independence test.

Scores used:
- Factorial Medium Original, Stability-Selected, Age-Only, and PC1 from `outputs_final/facs_v2_1/model_comparison/model_comparison_nested_scores.csv`.
- Elastic Net from `outputs_final/facs_v2_1/baselines/elastic_net/facs_v2_1_elastic_net_heldout_score_points.csv`, using held-out predicted Young probability as the youth-oriented score.

Metrics:
- AUC: Young vs Old, higher score indicates more Young-like.
- all-age rho: Spearman correlation between score and age months across all mice.
- old-only rho: Spearman correlation within 18m/24m Old mice only.
- library rho: Spearman correlation between score and pseudobulk library size.
- young-old median difference: median(Young score) - median(Old score).

Because the cohort has only 14 donors, bootstrap CIs should be interpreted as a small-sample uncertainty audit, not a substitute for external validation.

## Bootstrap CI Summary

| Model | Metric | Observed | Bootstrap median | 95% CI | Valid resamples |
|---|---:|---:|---:|---:|---:|
| factorial_medium_original | auc | 0.979 | 1.000 | [0.878, 1.000] | 9999/10000 |
| factorial_medium_original | all_age_rho | -0.633 | -0.654 | [-0.838, -0.120] | 10000/10000 |
| factorial_medium_original | old_only_rho | 0.655 | 0.677 | [0.000, 0.910] | 9815/10000 |
| factorial_medium_original | library_rho | -0.521 | -0.523 | [-0.794, 0.100] | 10000/10000 |
| factorial_medium_original | young_old_median_difference | 0.263 | 0.260 | [0.155, 0.395] | 9999/10000 |
| factorial_stability_selected | auc | 1.000 | 1.000 | [1.000, 1.000] | 9990/10000 |
| factorial_stability_selected | all_age_rho | -0.806 | -0.810 | [-0.923, -0.465] | 10000/10000 |
| factorial_stability_selected | old_only_rho | -0.109 | -0.120 | [-0.877, 0.828] | 9841/10000 |
| factorial_stability_selected | library_rho | -0.692 | -0.679 | [-0.920, -0.219] | 10000/10000 |
| factorial_stability_selected | young_old_median_difference | 0.260 | 0.256 | [0.108, 0.443] | 9990/10000 |
| age_only_de | auc | 1.000 | 1.000 | [1.000, 1.000] | 9997/10000 |
| age_only_de | all_age_rho | -0.769 | -0.780 | [-0.904, -0.364] | 10000/10000 |
| age_only_de | old_only_rho | 0.109 | 0.110 | [-0.784, 0.900] | 9813/10000 |
| age_only_de | library_rho | -0.727 | -0.705 | [-0.938, -0.277] | 10000/10000 |
| age_only_de | young_old_median_difference | 0.264 | 0.265 | [0.151, 0.413] | 9997/10000 |
| pc1_baseline | auc | 0.854 | 0.875 | [0.583, 1.000] | 9998/10000 |
| pc1_baseline | all_age_rho | -0.511 | -0.520 | [-0.829, 0.026] | 10000/10000 |
| pc1_baseline | old_only_rho | 0.218 | 0.234 | [-0.707, 0.900] | 9820/10000 |
| pc1_baseline | library_rho | -0.424 | -0.412 | [-0.774, 0.129] | 10000/10000 |
| pc1_baseline | young_old_median_difference | 0.221 | 0.221 | [-0.004, 0.525] | 9998/10000 |
| elastic_net | auc | 0.646 | 0.667 | [0.286, 1.000] | 9996/10000 |
| elastic_net | all_age_rho | -0.206 | -0.211 | [-0.884, 0.408] | 10000/10000 |
| elastic_net | old_only_rho | -0.109 | -0.133 | [-0.880, 0.811] | 9838/10000 |
| elastic_net | library_rho | -0.366 | -0.381 | [-0.864, 0.261] | 10000/10000 |
| elastic_net | young_old_median_difference | 0.054 | 0.054 | [-0.077, 0.217] | 9996/10000 |

## Interpretation Guardrails

- Wide or boundary-touching intervals are expected with 14 donors.
- Old-only rho has fewer informative resamples because resampled Old donors may lack both 18m and 24m age values.
- Library rho is reported as a diagnostic technical association; it is not corrected by this bootstrap.
- These intervals summarize the already frozen nested-score outputs and do not rerun feature selection inside each bootstrap resample.
