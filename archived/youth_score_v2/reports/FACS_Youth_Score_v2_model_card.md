# FACS Limb Muscle MSC Youth Score v2 Model Card

## Model Role

Primary frozen model: `factorial_medium_original` with 50 young-high and 50 old-high genes. Comparators exported in the same artifact set: `factorial_large_original`, `factorial_stability_selected`, `factorial_medium_equal_weight`, and `age_only_de`.

This is a FACS TMS Limb Muscle MSC cohort-state score. It should not be described as an externally validated universal MSC aging clock.

## Training Data

- Donors: 14 mice
- Cells represented by pseudobulk: 935
- Ages: 18m=4; 24m=4; 3m=6
- Sex by age group: =2; =2; =4; =6
- Filtered genes: 12859
- Full factorial design full-rank: TRUE

## Frozen Training Procedure

1. Filter genes at CPM > 1 in at least 2 mice.
2. Use edgeR TMM and limma voom for DE/ranking.
3. Fit full factorial age-by-sex group contrasts and age-only DE.
4. For the primary model, require female/male direction concordance, age-only/common direction concordance, exclusion of listed sex-linked genes, and finite common effects.
5. Rank by `abs(common_logFC) * abs(common_t) * interaction_penalty`.
6. Select top 50 young-high and top 50 old-high genes for Medium Original.
7. Freeze per-gene training mean and sd on full-data log2(CPM + 1).
8. Calibrate raw module difference with full training medians: `(raw - old_reference_center) / (young_reference_center - old_reference_center)`.

Training-set Young and Old medians are therefore fixed to 1 and 0 by construction; apparent full-data separation is not validation evidence.

## Validation Evidence Already Completed

- Donor bootstrap quantified uncertainty for primary and baseline models.
- Formal 999 sex-stratified complete-age-label nested permutation was run for the frozen primary pipeline.
- The primary permutation statistic was `abs_all_age_rho`; AUC distance from 0.5 was supporting only.

## Key Final Internal Results

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

## Donor Bootstrap: Primary Model

                     model                      metric   observed boot_median
 factorial_medium_original                         auc  1.0000000   1.0000000
 factorial_medium_original                 all_age_rho -0.7125870  -0.7313103
 factorial_medium_original                old_only_rho  0.4364358   0.4642383
 factorial_medium_original                 library_rho -0.8461538  -0.8318584
 factorial_medium_original young_old_median_difference  0.2648294   0.2648294
     ci_low    ci_high valid_bootstrap_n total_bootstrap_n
  1.0000000  1.0000000              9999             10000
 -0.8702076 -0.2383778             10000             10000
 -0.4277926  0.9065968              9815             10000
 -0.9864253 -0.5590776             10000             10000
  0.1499612  0.4385786              9999             10000

## Formal Permutation: Primary Model

                 metric                 role  observed null_median null_ci_low
        abs_all_age_rho              primary 0.7125870  0.33754121  0.01406422
 auc_distance_from_half           supporting 0.5000000  0.20833333  0.00000000
 young_minus_old_median recorded_not_primary 0.2648294 -0.07482546 -0.29566505
 null_ci_high empirical_p                         p_mode valid_null_n
    0.8063484       0.103        upper_tail_abs_spearman          999
    0.5000000       0.036   upper_tail_abs_auc_minus_0.5          999
    0.1662452       0.006 upper_tail_positive_difference          999

## Parser Validation

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

## Limitations

- Strong score-level library-size association remains in the FACS cohort.
- The 3m/18m frozen extension to 24m did not support a monotonic continuous aging trajectory.
- The cohort has only 14 mice and unbalanced 24m sex composition.
- Nested permutation tests label association under the predefined training pipeline, not biological specificity or technical independence.
- No independent external dataset validation has been completed yet.

## Exported Files

- `models/facs_v2_full_data_frozen_signatures_all_models.csv`
- `models/facs_v2_full_data_frozen_calibration.csv`
- `models/facs_v2_full_data_frozen_models.rds`
- `parser/score_facs_v2_youth_model.R`
- `full_data_frozen_training_scores.csv`
- `parser_numerical_equivalence.csv`
