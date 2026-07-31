# FACS v2 Step 14: Main and Comparator Model Runs

## Scope

Technical-tightening is closed. This step runs predefined main/comparator models without retuning thresholds, signature sizes, penalty strengths, or donor exclusions.

Models run here:

- Factorial Medium Original
- Factorial Large Original
- Factorial Medium Equal Weight
- Factorial Stability-Selected using inner-fold >=75% selection inside each outer fold
- Age-Only Signed-DE baseline
- Nested PC1 baseline

Kaile frozen model was not run because no Kaile model file is currently available. GSE176206 was not used.

No files under `data_facs` were modified.

## Summary

                         model       auc all_age_rho old_only_rho
     factorial_medium_original 1.0000000  -0.7125870    0.4364358
      factorial_large_original 1.0000000  -0.6938347    0.5455447
 factorial_medium_equal_weight 1.0000000  -0.6750824    0.6546537
  factorial_stability_selected 1.0000000  -0.7688439    0.1091089
                   age_only_de 1.0000000  -0.7500916    0.2182179
                  pc1_baseline 0.8333333  -0.4125504    0.3273268
 young_minus_old_median library_rho detected_rho cell_count_rho
              0.2648294  -0.8461538  -0.04175824     -0.7106715
              0.2195271  -0.8417582  -0.03296703     -0.6754680
              0.2442443  -0.8461538  -0.01098901     -0.6380642
              0.2821457  -0.8329670  -0.12087912     -0.7744779
              0.3150379  -0.8813187   0.01538462     -0.6732677
              0.2206597  -0.3450549   0.38021978      0.0924093
 median_signature_size median_jaccard genes_selected_ge_75pct
                 100.0      0.4492754                      57
                 200.0      0.4545455                     123
                 100.0      0.4492754                      57
                  56.5      0.4400000                      32
                 100.0      0.4705882                      52
                    NA             NA                      NA
 zero_signature_folds                    role
                    0      Primary predefined
                    0        Large comparator
                    0      Weight sensitivity
                    0    Stability comparator
                    0     Biological baseline
                    0 Technical-axis baseline

## Notes

- The PC1 model is a technical-axis comparator, not a biological Youth Score.
- The stability-selected comparator uses inner stability within each outer fold to avoid held-out leakage.
- Model roles remain predefined; this table is for credibility assessment, not post-hoc primary-model replacement.

## Outputs

- `model_comparison_nested_scores.csv`
- `model_comparison_fold_signatures.csv`
- `model_comparison_fold_diagnostics.csv`
- `model_comparison_summary.csv`
- `elastic_net_reference_summary.csv`
