# FACS v2.1 Step 13: Equal-Cell Subsampling Feasibility Audit

## Scope

This step follows the negative technical-penalty result. It asks whether donor-level score dependence on library size/cell count persists when each mouse contributes the same number of cells.

- Replicates per scenario: 20
- Full retraining mode: every replicate reruns fully nested LOMO Medium signature.
- Frozen rescoring mode: uses a frozen full-data Medium signature and rescales equal-cell pseudobulk.
- No files under `data_facs` were modified.

## Scenarios

- `all_mice_14_cells`: all 14 mice, 14 cells per mouse.
- `min25_mice_28_cells`: excludes mice with <25 cells, then uses 28 cells per mouse.

## Summary

            scenario                            mode n_replicates auc_median
   all_mice_14_cells frozen_full_data_medium_rescore           20  1.0000000
   all_mice_14_cells     full_retraining_nested_lomo           20  0.9479167
 min25_mice_28_cells frozen_full_data_medium_rescore           20  1.0000000
 min25_mice_28_cells     full_retraining_nested_lomo           20  0.7500000
  auc_q025  auc_q975 all_age_rho_median all_age_rho_q025 all_age_rho_q975
 1.0000000 1.0000000         -0.7875961       -0.9019851       -0.7304017
 0.7916667 1.0000000         -0.7430595       -0.8800684       -0.4200513
 1.0000000 1.0000000         -0.7982717       -0.9460998       -0.5868776
 0.4820313 0.9851562         -0.3400046       -0.7842281        0.1382193
 old_only_rho_median old_only_rho_q025 old_only_rho_q975
           0.0000000        -0.6655646         0.3327823
          -0.1636634        -0.6028269         0.6028269
          -0.3273268        -0.8728716         0.4528021
          -0.1091089        -0.4364358         0.7173913
 young_minus_old_median_median young_minus_old_median_q025
                     1.1007444                  0.91121357
                     0.1548642                  0.08060162
                     0.9934025                  0.84266146
                     0.1135995                 -0.01291111
 young_minus_old_median_q975 library_rho_median library_rho_q025
                   1.2776497         -0.5296703       -0.7676923
                   0.2821933         -0.4373626       -0.6312088
                   1.1199403         -0.5769231       -0.7853147
                   0.2460105         -0.3496503       -0.6017483
 library_rho_q975 detected_rho_median detected_rho_q025 detected_rho_q975
      -0.06021978           0.5956044        0.52945055         0.7695604
      -0.12527473           0.6329670        0.25648352         0.8356044
      -0.23339161           0.6118881        0.33269231         0.7669580
       0.17220280           0.5419580       -0.02954545         0.6425724
 cell_count_rho_median cell_count_rho_q025 cell_count_rho_q975
                    NA                  NA                  NA
                    NA                  NA                  NA
                    NA                  NA                  NA
                    NA                  NA                  NA
 zero_signature_folds_median zero_signature_folds_q025
                           0                         0
                           0                         0
                           0                         0
                           0                         0
 zero_signature_folds_q975 prop_auc_ge_0_8 prop_all_age_rho_negative
                         0             1.0                       1.0
                         0             0.9                       1.0
                         0             1.0                       1.0
                         0             0.4                       0.9
 prop_old_only_rho_negative
                       0.45
                       0.60
                       0.75
                       0.65

## Interpretation Guardrails

This is a 20-replicate feasibility run, not the final 100-500 replicate uncertainty analysis.
If the summary shows persistent high technical correlations or unstable age direction, the cohort-level technical structure remains unresolved.

## Outputs

- `equal_cell_subsampling_scores.csv`
- `equal_cell_subsampling_nested_signatures.csv`
- `equal_cell_subsampling_replicate_metrics.csv`
- `equal_cell_subsampling_summary.csv`
