# FACS v2.1 Step 12: Technical-Burden Sensitivity Audit

## Scope

The original factorial ranking is retained as the primary specification. Technical-burden penalties are sensitivity branches only.

Technical burden is estimated within each outer training fold from residual expression after `~0 + sex:age_group` adjustment. The held-out mouse is not used to estimate burden, ranking, weights, or calibration.

No files under `data_facs` were modified.

## Penalty

Composite burden: mean absolute residual Spearman correlation with raw library size, cell count, and detected genes.

Penalty: `1 / (1 + lambda * burden)`.

Branches: original lambda=0, tech_mild lambda=1, tech_moderate lambda=2, tech_strong lambda=4.

Lowest-depth donor full rerun excludes: `3_9_M`.

## Summary

                  analysis ranking_branch lambda signature_version n_folds
                  all_mice       original      0             Small      14
                  all_mice       original      0            Medium      14
                  all_mice       original      0             Large      14
                  all_mice      tech_mild      1             Small      14
                  all_mice      tech_mild      1            Medium      14
                  all_mice      tech_mild      1             Large      14
                  all_mice  tech_moderate      2             Small      14
                  all_mice  tech_moderate      2            Medium      14
                  all_mice  tech_moderate      2             Large      14
                  all_mice    tech_strong      4             Small      14
                  all_mice    tech_strong      4            Medium      14
                  all_mice    tech_strong      4             Large      14
 remove_lowest_depth_3_9_M       original      0             Small      13
 remove_lowest_depth_3_9_M       original      0            Medium      13
 remove_lowest_depth_3_9_M       original      0             Large      13
 remove_lowest_depth_3_9_M      tech_mild      1             Small      13
 remove_lowest_depth_3_9_M      tech_mild      1            Medium      13
 remove_lowest_depth_3_9_M      tech_mild      1             Large      13
 remove_lowest_depth_3_9_M  tech_moderate      2             Small      13
 remove_lowest_depth_3_9_M  tech_moderate      2            Medium      13
 remove_lowest_depth_3_9_M  tech_moderate      2             Large      13
 remove_lowest_depth_3_9_M    tech_strong      4             Small      13
 remove_lowest_depth_3_9_M    tech_strong      4            Medium      13
 remove_lowest_depth_3_9_M    tech_strong      4             Large      13
 zero_signature_folds       auc all_age_rho old_only_rho young_minus_old_median
                    0 1.0000000  -0.7500916    0.2182179              0.3137950
                    0 0.9791667  -0.6328898    0.6546537              0.2631161
                    0 1.0000000  -0.7125870    0.4364358              0.2827920
                    0 0.9791667  -0.6891466    0.3273268              0.3135429
                    0 0.9791667  -0.6328898    0.6546537              0.2426288
                    0 1.0000000  -0.7125870    0.4364358              0.2906874
                    0 0.9375000  -0.5860090    0.4364358              0.2908689
                    0 1.0000000  -0.6750824    0.6546537              0.2342594
                    0 1.0000000  -0.6938347    0.5455447              0.2797898
                    0 0.8958333  -0.4828714    0.5455447              0.3126447
                    0 1.0000000  -0.6750824    0.6546537              0.2481607
                    0 1.0000000  -0.7125870    0.4364358              0.2842534
                    0 0.7500000  -0.3319569    0.1091089              0.1863894
                    0 0.9000000  -0.5124948    0.3273268              0.2516988
                    0 1.0000000  -0.6639137    0.4364358              0.2728983
                    0 0.8500000  -0.4134901    0.4364358              0.2609750
                    0 0.9750000  -0.5911162    0.5455447              0.2470378
                    0 1.0000000  -0.6639137    0.4364358              0.2727877
                    0 0.8750000  -0.5095829    0.2182179              0.2292477
                    0 0.9750000  -0.5911162    0.5455447              0.2245506
                    0 1.0000000  -0.6406185    0.5455447              0.2433873
                    0 0.8500000  -0.3901949    0.5455447              0.2235758
                    0 0.9750000  -0.6377066    0.3273268              0.2609135
                    0 1.0000000  -0.6639137    0.4364358              0.2546232
 library_rho detected_rho cell_count_rho median_jaccard genes_selected_ge_75pct
  -0.6527473  -0.03736264     -0.5544558      0.3793103                      19
  -0.5208791  -0.11648352     -0.4620465      0.4285714                      58
  -0.6439560   0.17362637     -0.4070410      0.4388489                     122
  -0.6263736  -0.08571429     -0.5918595      0.4035088                      20
  -0.5604396  -0.01538462     -0.4488452      0.4184397                      55
  -0.6747253   0.23956044     -0.3806383      0.4440433                     127
  -0.5692308  -0.08131868     -0.4884491      0.3793103                      19
  -0.6087912   0.06813187     -0.3960398      0.4184397                      57
  -0.6791209   0.24395604     -0.3322334      0.4492754                     124
  -0.5032967  -0.03296703     -0.4136416      0.3559322                      18
  -0.5956044   0.08131868     -0.3894392      0.3986014                      53
  -0.6835165   0.23076923     -0.3520354      0.4440433                     121
  -0.3351648  -0.35164835     -0.7730406      0.3559322                      18
  -0.4945055  -0.23076923     -0.6437420      0.3937451                      54
  -0.6318681  -0.01648352     -0.4896841      0.4440433                     126
  -0.4230769  -0.27472527     -0.5612110      0.3559322                      19
  -0.5549451  -0.08241758     -0.4869331      0.3937451                      58
  -0.6373626   0.04945055     -0.4401655      0.4492754                     124
  -0.4835165  -0.29120879     -0.7372772      0.3446328                      18
  -0.5384615  -0.03296703     -0.4209082      0.4084507                      57
  -0.6263736   0.07692308     -0.3961489      0.4492754                     124
  -0.3626374  -0.25274725     -0.5694641      0.3559322                      19
  -0.5659341  -0.06593407     -0.5281986      0.3793760                      53
  -0.6758242   0.06593407     -0.4731779      0.4285714                     115
 median_signature_technical_burden
                         0.2067258
                         0.1985396
                         0.2003086
                         0.1655983
                         0.1722842
                         0.1770990
                         0.1488638
                         0.1574421
                         0.1689646
                         0.1336600
                         0.1443290
                         0.1574285
                         0.2078726
                         0.2065199
                         0.2121212
                         0.1789511
                         0.1946409
                         0.2009754
                         0.1730697
                         0.1827963
                         0.1939164
                         0.1577259
                         0.1717490
                         0.1772277

## Outputs

- `technical_burden_gene_metrics_by_fold.csv`
- `technical_penalty_fold_signatures.csv`
- `technical_penalty_nested_lomo_scores.csv`
- `technical_penalty_fold_diagnostics.csv`
- `technical_penalty_model_summary.csv`
- `technical_penalty_gene_selection_frequency.csv`
