# FACS v2 Step 13: Equal-Cell Subsampling Feasibility Audit

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
   all_mice_14_cells frozen_full_data_medium_rescore           20    1.00000
   all_mice_14_cells     full_retraining_nested_lomo           20    1.00000
 min25_mice_28_cells frozen_full_data_medium_rescore           20    1.00000
 min25_mice_28_cells     full_retraining_nested_lomo           20    0.90625
  auc_q025 auc_q975 all_age_rho_median all_age_rho_q025 all_age_rho_q975
 0.9890625        1         -0.7875961       -0.8724503       -0.7190331
 0.9166667        1         -0.7407154       -0.8654181       -0.5892907
 1.0000000        1         -0.7243577       -0.8744032       -0.6193997
 0.7648438        1         -0.5321812       -0.8152720       -0.2365250
 old_only_rho_median old_only_rho_q025 old_only_rho_q975
          0.00000000        -0.4937180         0.3273268
          0.05455447        -0.5564556         0.4364358
         -0.05455447        -0.6082824         0.3327823
          0.16366342        -0.4418912         0.6546537
 young_minus_old_median_median young_minus_old_median_q025
                     1.1413544                  0.98368721
                     0.2624977                  0.15439120
                     1.1068948                  0.97994180
                     0.1845504                  0.06425355
 young_minus_old_median_q975 library_rho_median library_rho_q025
                   1.3173498         -0.5736264       -0.7979121
                   0.3376360         -0.6021978       -0.7371429
                   1.1881955         -0.6608392       -0.8648601
                   0.2935824         -0.5944056       -0.8365385
 library_rho_q975 detected_rho_median detected_rho_q025 detected_rho_q975
       -0.2892308           0.6109890         0.4430769         0.8125275
       -0.3215385           0.6241758         0.3538462         0.7762637
       -0.5237762           0.5000000         0.3173077         0.7017483
       -0.2211538           0.3321678         0.1291958         0.5701049
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
                         0             1.0                         1
                         0             1.0                         1
                         0             1.0                         1
                         0             0.9                         1
 prop_old_only_rho_negative
                       0.45
                       0.35
                       0.50
                       0.30

## Interpretation Guardrails

This is a 20-replicate feasibility run, not the final 100-500 replicate uncertainty analysis.
If the summary shows persistent high technical correlations or unstable age direction, the cohort-level technical structure remains unresolved.

## Outputs

- `equal_cell_subsampling_scores.csv`
- `equal_cell_subsampling_nested_signatures.csv`
- `equal_cell_subsampling_replicate_metrics.csv`
- `equal_cell_subsampling_summary.csv`
