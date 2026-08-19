# FACS v2.1 Step 14: Main and Comparator Model Runs

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
     factorial_medium_original 0.9791667  -0.6328898    0.6546537
      factorial_large_original 1.0000000  -0.7125870    0.4364358
 factorial_medium_equal_weight 0.9791667  -0.6141375    0.7637626
  factorial_stability_selected 1.0000000  -0.8063484   -0.1091089
                   age_only_de 1.0000000  -0.7688439    0.1091089
                  pc1_baseline 0.8541667  -0.5109999    0.2182179
 young_minus_old_median library_rho detected_rho cell_count_rho
              0.2631161  -0.5208791  -0.11648352    -0.46204648
              0.2827920  -0.6439560   0.17362637    -0.40704095
              0.2576940  -0.5428571  -0.01978022    -0.37843807
              0.2601991  -0.6923077   0.04175824    -0.55665600
              0.2637327  -0.7274725   0.13406593    -0.51265158
              0.2214713  -0.4241758   0.38901099     0.02860288
 median_signature_size median_jaccard genes_selected_ge_75pct
                 100.0      0.4285714                      58
                 200.0      0.4388489                     122
                 100.0      0.4285714                      58
                  54.5      0.4074074                      30
                 100.0      0.4492754                      54
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
