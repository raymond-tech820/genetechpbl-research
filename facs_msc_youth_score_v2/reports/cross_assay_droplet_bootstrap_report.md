# Cross-Assay Droplet Donor Bootstrap Uncertainty

Bootstrap iterations: 10000 donor-level resamples with replacement per model.

Purpose: quantify uncertainty for the already frozen cross-assay sensitivity results. This does not retrain models, tune thresholds, or select a new primary model.

Primary focus:
- `FACS_v2__factorial_medium_original`: FACS-trained frozen primary applied to Droplet.
- `Droplet_v1__Medium`: Droplet-trained frozen primary applied to Droplet as same-assay reference.

## Primary CI Summary

| Model | Metric | Observed | Bootstrap median | 95% CI | Valid resamples |
|---|---:|---:|---:|---:|---:|
| Droplet_v1__Medium | auc_young_vs_old | 1.000 | 1.000 | [1.000, 1.000] | 8824/10000 |
| Droplet_v1__Medium | all_age_rho | -0.939 | -0.933 | [-0.979, -0.774] | 10000/10000 |
| Droplet_v1__Medium | old_only_rho | -0.895 | -0.885 | [-0.964, -0.639] | 9998/10000 |
| Droplet_v1__Medium | raw_library_rho | 0.133 | 0.143 | [-0.445, 0.619] | 10000/10000 |
| Droplet_v1__Medium | effective_library_rho | 0.098 | 0.107 | [-0.489, 0.594] | 10000/10000 |
| Droplet_v1__Medium | detected_genes_rho | 0.182 | 0.180 | [-0.458, 0.714] | 10000/10000 |
| Droplet_v1__Medium | cell_count_rho | 0.063 | 0.071 | [-0.491, 0.585] | 10000/10000 |
| Droplet_v1__Medium | young_minus_old_median | 1.000 | 1.000 | [0.762, 1.057] | 8824/10000 |
| FACS_v2__factorial_medium_original | auc_young_vs_old | 1.000 | 1.000 | [1.000, 1.000] | 8873/10000 |
| FACS_v2__factorial_medium_original | all_age_rho | -0.873 | -0.881 | [-0.971, -0.538] | 10000/10000 |
| FACS_v2__factorial_medium_original | old_only_rho | -0.778 | -0.794 | [-0.951, -0.233] | 9993/10000 |
| FACS_v2__factorial_medium_original | raw_library_rho | -0.238 | -0.239 | [-0.707, 0.400] | 10000/10000 |
| FACS_v2__factorial_medium_original | effective_library_rho | -0.259 | -0.262 | [-0.736, 0.372] | 10000/10000 |
| FACS_v2__factorial_medium_original | detected_genes_rho | 0.119 | 0.120 | [-0.548, 0.696] | 10000/10000 |
| FACS_v2__factorial_medium_original | cell_count_rho | -0.308 | -0.307 | [-0.717, 0.314] | 10000/10000 |
| FACS_v2__factorial_medium_original | young_minus_old_median | 0.150 | 0.150 | [0.062, 0.226] | 8873/10000 |

## Interpretation Guardrails

- The Droplet cohort has only 12 mice, including 2 Young donors; wide and boundary-touching CIs are expected.
- AUC bootstrap resamples may be invalid when a resample lacks Young or Old donors; valid resample counts are reported.
- Old-only rho is based on 18m/21m/24m Old mice and remains a small-sample diagnostic.
- Technical-variable CIs quantify association in the Droplet application; they do not prove technical independence.

## Outputs

- `cross_assay_droplet_bootstrap_ci.csv`
- `cross_assay_droplet_primary_bootstrap_ci.csv`
- `cross_assay_droplet_bootstrap_samples.csv`
- `cross_assay_droplet_bootstrap_observed_metrics.csv`
