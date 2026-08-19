# FACS Youth Score v2.1 Step 18: Full-Data Frozen Model Export

## Scope

This step freezes the finalized FACS v2.1 model artifacts. It does not retune thresholds, change model roles, or inspect cross-assay performance.

Primary: `factorial_medium_original`.

Comparators: `factorial_large_original`, `factorial_stability_selected`, `factorial_medium_equal_weight`, `age_only_de`.

## Model Sizes

                         model signature_size young_high_n old_high_n
     factorial_medium_original            100           50         50
      factorial_large_original            200          100        100
  factorial_stability_selected             58           30         28
 factorial_medium_equal_weight            100           50         50
                   age_only_de            100           50         50
 calibration_denominator
                2.545538
                2.458266
                2.650305
                2.510409
                2.738931

## Full-Data Apparent Scores

                         model auc_young_vs_old all_age_rho old_only_rho
                   age_only_de                1  -0.8251007   -0.2182179
      factorial_large_original                1  -0.7313393    0.3273268
 factorial_medium_equal_weight                1  -0.7313393    0.3273268
     factorial_medium_original                1  -0.7313393    0.3273268
  factorial_stability_selected                1  -0.7688439    0.1091089
 young_minus_old_median library_rho detected_rho cell_count_rho
                      1  -0.7890110   0.27032967     -0.4906494
                      1  -0.6747253  -0.01098901     -0.5258529
                      1  -0.6659341   0.04175824     -0.5258529
                      1  -0.6659341   0.04175824     -0.5258529
                      1  -0.7098901  -0.09450549     -0.6908695

These are apparent full-training scores. They are useful for artifact verification and calibration, not validation evidence.

## Parser Numerical Equivalence

                         model n_samples min_gene_coverage
     factorial_medium_original        14                 1
      factorial_large_original        14                 1
  factorial_stability_selected        14                 1
 factorial_medium_equal_weight        14                 1
                   age_only_de        14                 1
 min_weighted_coverage max_abs_raw_score_diff max_abs_calibrated_score_diff
                     1           1.332268e-15                  1.050375e-15
                     1           1.332268e-15                  3.330669e-15
                     1           1.332268e-15                  1.110223e-15
                     1           1.110223e-15                  2.442491e-15
                     1           2.664535e-15                  1.110223e-15

## Stability-Selected Definition

Stable pool threshold: selected in at least 75% of Step 14 outer folds. Stable pool size: 58 genes. Exported stability-selected signature size: 58 genes.

## Outputs

- `models/facs_v2_1_full_data_frozen_signatures_all_models.csv`
- `models/facs_v2_1_full_data_frozen_calibration.csv`
- `models/facs_v2_1_full_data_frozen_models.rds`
- `parser/score_facs_v2_1_youth_model.R`
- `full_data_frozen_training_scores.csv`
- `full_data_frozen_training_score_metrics.csv`
- `parser_numerical_equivalence.csv`
- `FACS_Youth_Score_v2_1_model_card.md`
