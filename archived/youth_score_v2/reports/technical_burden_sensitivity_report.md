# FACS v2 Step 12: Technical-Burden Sensitivity Audit

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
 zero_signature_folds   auc all_age_rho old_only_rho young_minus_old_median
                    0 1.000  -0.7688439    0.1091089              0.3362587
                    0 1.000  -0.7125870    0.4364358              0.2648294
                    0 1.000  -0.6938347    0.5455447              0.2195271
                    0 1.000  -0.7688439    0.1091089              0.3008770
                    0 1.000  -0.7313393    0.3273268              0.2205958
                    0 1.000  -0.6938347    0.5455447              0.2212418
                    0 1.000  -0.7875961    0.0000000              0.3644530
                    0 1.000  -0.7313393    0.3273268              0.2378911
                    0 1.000  -0.7125870    0.4364358              0.2672098
                    0 1.000  -0.7500916    0.2182179              0.3076009
                    0 1.000  -0.7688439    0.1091089              0.2372376
                    0 1.000  -0.6938347    0.5455447              0.2746151
                    0 0.900  -0.5357900    0.3273268              0.1834742
                    0 0.950  -0.6114995    0.3273268              0.2416141
                    0 1.000  -0.6872089    0.3273268              0.1874410
                    0 0.925  -0.6085876    0.3273268              0.1324559
                    0 0.950  -0.6114995    0.3273268              0.2547938
                    0 1.000  -0.6639137    0.4364358              0.1626964
                    0 0.925  -0.5852924    0.3273268              0.2346004
                    0 0.950  -0.6347947    0.2182179              0.2536383
                    0 1.000  -0.6639137    0.4364358              0.1931535
                    0 0.975  -0.6610018    0.3273268              0.1720906
                    0 0.950  -0.6347947    0.2182179              0.1846652
                    0 1.000  -0.6639137    0.4364358              0.1804348
 library_rho detected_rho cell_count_rho median_jaccard genes_selected_ge_75pct
  -0.8241758  -0.12967033     -0.7678773      0.4545455                      24
  -0.8461538  -0.04175824     -0.7106715      0.4492754                      57
  -0.8417582  -0.03296703     -0.6754680      0.4545455                     123
  -0.8241758  -0.14285714     -0.7502755      0.4545455                      24
  -0.8681319   0.05494505     -0.6578662      0.4285714                      57
  -0.8593407  -0.01098901     -0.6600664      0.4545455                     117
  -0.8373626  -0.09890110     -0.7832788      0.4545455                      25
  -0.8637363   0.03296703     -0.6622666      0.4285714                      54
  -0.8461538   0.02417582     -0.6490653      0.4545455                     121
  -0.8109890  -0.13406593     -0.7216726      0.4285714                      22
  -0.8769231   0.05494505     -0.7172722      0.4084507                      51
  -0.8769231  -0.04615385     -0.6292633      0.4336918                     112
  -0.7472527  -0.23626374     -0.7345261      0.3793103                      19
  -0.7967033  -0.12637363     -0.7647875      0.4035261                      57
  -0.8516484  -0.06593407     -0.7235220      0.4234875                     123
  -0.7307692  -0.18681319     -0.7317751      0.3676213                      21
  -0.8021978  -0.12087912     -0.7565344      0.4085206                      61
  -0.8626374  -0.01648352     -0.6547461      0.4235056                     118
  -0.7472527  -0.25824176     -0.7812937      0.3793103                      19
  -0.8241758  -0.09340659     -0.7757917      0.4084507                      58
  -0.8626374  -0.01648352     -0.6547461      0.4285714                     120
  -0.8241758  -0.17582418     -0.7317751      0.3793103                      22
  -0.8296703  -0.08791209     -0.7675386      0.3888889                      55
  -0.8571429  -0.02197802     -0.6629993      0.4159336                     104
 median_signature_technical_burden
                         0.1785013
                         0.1833202
                         0.1848940
                         0.1509982
                         0.1622453
                         0.1699581
                         0.1398631
                         0.1462249
                         0.1600736
                         0.1260484
                         0.1303056
                         0.1446374
                         0.1914863
                         0.1924466
                         0.1907259
                         0.1710341
                         0.1777944
                         0.1781907
                         0.1479582
                         0.1620993
                         0.1694174
                         0.1323794
                         0.1466500
                         0.1546210

## Outputs

- `technical_burden_gene_metrics_by_fold.csv`
- `technical_penalty_fold_signatures.csv`
- `technical_penalty_nested_lomo_scores.csv`
- `technical_penalty_fold_diagnostics.csv`
- `technical_penalty_model_summary.csv`
- `technical_penalty_gene_selection_frequency.csv`
