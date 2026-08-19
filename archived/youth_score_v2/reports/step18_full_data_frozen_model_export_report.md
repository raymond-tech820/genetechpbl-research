# FACS Youth Score v2 Step 18: Full-Data Frozen Model Export

## Scope

This step freezes the finalized FACS v2 model artifacts. It does not retune thresholds, change model roles, or inspect cross-assay performance.

Primary: `factorial_medium_original`.

Comparators: `factorial_large_original`, `factorial_stability_selected`, `factorial_medium_equal_weight`, `age_only_de`.

## Model Sizes

                         model signature_size young_high_n old_high_n
     factorial_medium_original            100           50         50
      factorial_large_original            200          100        100
  factorial_stability_selected             57           29         28
 factorial_medium_equal_weight            100           50         50
                   age_only_de            100           50         50
 calibration_denominator
                2.388366
                2.351519
                2.558554
                2.320877
                2.720908

## Full-Data Apparent Scores

                         model auc_young_vs_old all_age_rho old_only_rho
                   age_only_de                1  -0.7688439    0.1091089
      factorial_large_original                1  -0.7313393    0.3273268
 factorial_medium_equal_weight                1  -0.7500916    0.2182179
     factorial_medium_original                1  -0.7500916    0.2182179
  factorial_stability_selected                1  -0.7875961    0.0000000
 young_minus_old_median library_rho detected_rho cell_count_rho
                      1  -0.8505495  -0.05934066     -0.7480753
                      1  -0.8153846  -0.01538462     -0.6600664
                      1  -0.7978022   0.01098901     -0.6534657
                      1  -0.8109890  -0.03296703     -0.6886693
                      1  -0.8109890   0.08571429     -0.6468651

These are apparent full-training scores. They are useful for artifact verification and calibration, not validation evidence.

## Parser Numerical Equivalence

                         model n_samples min_gene_coverage
     factorial_medium_original        14                 1
      factorial_large_original        14                 1
  factorial_stability_selected        14                 1
 factorial_medium_equal_weight        14                 1
                   age_only_de        14                 1
 min_weighted_coverage max_abs_raw_score_diff max_abs_calibrated_score_diff
                     1           1.554312e-15                  1.776357e-15
                     1           1.554312e-15                  1.332268e-15
                     1           2.220446e-15                  2.331468e-15
                     1           6.661338e-16                  6.661338e-16
                     1           2.220446e-15                  4.857226e-16

## Stability-Selected Definition

Stable pool threshold: selected in at least 75% of Step 14 outer folds. Stable pool size: 57 genes. Exported stability-selected signature size: 57 genes.

## Outputs

- `models/facs_v2_full_data_frozen_signatures_all_models.csv`
- `models/facs_v2_full_data_frozen_calibration.csv`
- `models/facs_v2_full_data_frozen_models.rds`
- `parser/score_facs_v2_youth_model.R`
- `full_data_frozen_training_scores.csv`
- `full_data_frozen_training_score_metrics.csv`
- `parser_numerical_equivalence.csv`
- `FACS_Youth_Score_v2_model_card.md`
