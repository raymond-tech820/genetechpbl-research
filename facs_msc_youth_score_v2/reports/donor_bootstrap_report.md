# FACS Youth Score v2 Step 16: Donor Bootstrap Uncertainty

Bootstrap iterations: 10000 donor-level resamples with replacement per model.

Purpose: quantify uncertainty for predefined candidate/baseline models. This analysis is not a winner-selection or technical-independence test.

Scores used:
- Factorial Medium Original, Stability-Selected, Age-Only, and PC1 from `outputs/facs_v2/model_comparison/model_comparison_nested_scores.csv`.
- Elastic Net from `outputs/facs_v2/baselines/elastic_net/facs_v2_elastic_net_heldout_score_points.csv`, using held-out predicted Young probability as the youth-oriented score.

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
| factorial_medium_original | auc | 1.000 | 1.000 | [1.000, 1.000] | 9999/10000 |
| factorial_medium_original | all_age_rho | -0.713 | -0.731 | [-0.870, -0.238] | 10000/10000 |
| factorial_medium_original | old_only_rho | 0.436 | 0.464 | [-0.428, 0.907] | 9815/10000 |
| factorial_medium_original | library_rho | -0.846 | -0.832 | [-0.986, -0.559] | 10000/10000 |
| factorial_medium_original | young_old_median_difference | 0.265 | 0.265 | [0.150, 0.439] | 9999/10000 |
| factorial_stability_selected | auc | 1.000 | 1.000 | [1.000, 1.000] | 9990/10000 |
| factorial_stability_selected | all_age_rho | -0.769 | -0.781 | [-0.901, -0.366] | 10000/10000 |
| factorial_stability_selected | old_only_rho | 0.109 | 0.129 | [-0.816, 0.877] | 9841/10000 |
| factorial_stability_selected | library_rho | -0.833 | -0.814 | [-0.977, -0.540] | 10000/10000 |
| factorial_stability_selected | young_old_median_difference | 0.282 | 0.306 | [0.166, 0.529] | 9990/10000 |
| age_only_de | auc | 1.000 | 1.000 | [1.000, 1.000] | 9997/10000 |
| age_only_de | all_age_rho | -0.750 | -0.768 | [-0.894, -0.318] | 10000/10000 |
| age_only_de | old_only_rho | 0.218 | 0.257 | [-0.725, 0.905] | 9813/10000 |
| age_only_de | library_rho | -0.881 | -0.869 | [-0.986, -0.599] | 10000/10000 |
| age_only_de | young_old_median_difference | 0.315 | 0.320 | [0.171, 0.489] | 9997/10000 |
| pc1_baseline | auc | 0.833 | 0.850 | [0.556, 1.000] | 9998/10000 |
| pc1_baseline | all_age_rho | -0.413 | -0.424 | [-0.847, 0.184] | 10000/10000 |
| pc1_baseline | old_only_rho | 0.327 | 0.329 | [-0.642, 0.905] | 9820/10000 |
| pc1_baseline | library_rho | -0.345 | -0.339 | [-0.673, 0.182] | 10000/10000 |
| pc1_baseline | young_old_median_difference | 0.221 | 0.221 | [-0.037, 0.536] | 9998/10000 |
| elastic_net | auc | 0.833 | 0.850 | [0.556, 1.000] | 9996/10000 |
| elastic_net | all_age_rho | -0.531 | -0.536 | [-0.873, -0.024] | 10000/10000 |
| elastic_net | old_only_rho | 0.000 | 0.000 | [-0.866, 0.854] | 9838/10000 |
| elastic_net | library_rho | -0.378 | -0.375 | [-0.756, 0.211] | 10000/10000 |
| elastic_net | young_old_median_difference | 0.224 | 0.247 | [-0.029, 0.461] | 9996/10000 |

## Interpretation Guardrails

- Wide or boundary-touching intervals are expected with 14 donors.
- Old-only rho has fewer informative resamples because resampled Old donors may lack both 18m and 24m age values.
- Library rho is reported as a diagnostic technical association; it is not corrected by this bootstrap.
- These intervals summarize the already frozen nested-score outputs and do not rerun feature selection inside each bootstrap resample.
