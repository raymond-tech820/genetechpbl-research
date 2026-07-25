# Droplet v1.1 vs v2.1 Formal Comparison Report

## Output Definitions And Aggregation

The two models produce young-like scores with the same direction but different physical meanings.

Kaile v1.1 scores each cell first. For cell $i$, the model outputs a calibrated probability-like score:

$$
p_i = P(Young \mid x_i)
$$

The donor score is the arithmetic mean of cell-level probabilities within mouse $m$:

$$
Y^{Kaile}_m = \frac{1}{n_m} \sum_{i \in m} p_i
$$

Thus Kaile v1.1 follows:

$$
single\ cells \rightarrow cell\ scores \rightarrow donor\ mean\ score
$$

Zihan v2.1 first aggregates raw single-cell counts into donor pseudobulk:

$$
C_{g,m} = \sum_{i \in m} c_{g,i}
$$

It then computes deployable expression and frozen gene standardization:

$$
E_{g,m} = \log_2\left(\frac{10^6 C_{g,m}}{\sum_h C_{h,m}} + 1\right)
$$

$$
Z_{g,m} = \frac{E_{g,m} - \mu_g}{s_g}
$$

Young-high and old-high module scores are weighted means:

$$
M^Y_m = \frac{\sum_{g \in G_Y} w_g Z_{g,m}}{\sum_{g \in G_Y} w_g}
$$

$$
M^O_m = \frac{\sum_{g \in G_O} w_g Z_{g,m}}{\sum_{g \in G_O} w_g}
$$

The raw and calibrated v2.1 donor scores are:

$$
S^{raw}_m = M^Y_m - M^O_m
$$

$$
Y^{Zihan}_m = \frac{S^{raw}_m - C_{old}}{C_{young} - C_{old}}
$$

Thus Zihan v2.1 follows:

$$
single\ cells \rightarrow donor\ pseudobulk \rightarrow donor\ score
$$

Both scores are oriented so that higher values indicate a more young-like transcriptional state. Their absolute numeric units are not equivalent, so comparisons should focus on donor ranking, age association, young-old separation, technical association, and concordance rather than equal score magnitudes.


## Scope

Both frozen FACS-derived models were applied to local `data_droplet/raw_data`. No model was retrained or modified. This is cross-assay sensitivity and transportability analysis, not external validation.

## Input Data

```text
https://drive.google.com/drive/folders/1wEx3eIoCxjNavELDaVkRrRuONPeqFiqp?usp=share_link
```

| Check | Value |
|---|---:|
| Cells | 9649 |
| Genes | 20138 |
| Mice | 12 |
| Diaphragm-like subtissue cells | 0 |
| Non-Limb_Muscle cells | 0 |
| Non-MSC cells | 0 |
| 1m cells | 0 |
| 30m cells | 0 |

This local Droplet raw_data contains 12 donors and excludes the 1m and 30m donors present in the packaged Kaile v1.1 Droplet sensitivity output.

## Packaged v1.1 Droplet Check

| Check | Value |
|---|---:|
| Overlap donors | 12 |
| Max absolute score difference on overlap | 5.210e-09 |
| Spearman local vs packaged on overlap | 1.000 |
| Local cells | 9649 |
| Packaged sensitivity cells | 13037 |
| Local donors | 12 |
| Packaged sensitivity donors | 16 |

## Frozen Cross-Assay Application

| Model | AUC | All-age Spearman | Old-only Spearman | Young-old median difference |
|---|---:|---:|---:|---:|
| Kaile v1.1 gene_signature | 0.950 | -0.873 | -0.817 | 0.126 |
| Zihan v2.1 factorial_medium_original | 1.000 | -0.830 | -0.701 | 0.125 |

## Paired Donor Bootstrap

| Metric | Kaile observed [95% bootstrap CI] | Zihan observed [95% bootstrap CI] | Delta Kaile minus Zihan [95% bootstrap CI] |
|---|---:|---:|---:|
| AUC | 0.950 [0.727, 1.000] | 1.000 [1.000, 1.000] | -0.050 [-0.273, 0.000] |
| All-age Spearman | -0.873 [-0.964, -0.624] | -0.830 [-0.954, -0.384] | -0.044 [-0.412, 0.208] |
| Old-only Spearman | -0.817 [-0.953, -0.427] | -0.701 [-0.895, -0.037] | -0.117 [-0.691, 0.342] |
| Young-old median difference | 0.126 [0.071, 0.176] | 0.125 [0.058, 0.174] | 0.002 [-0.033, 0.060] |

The Droplet point estimates show similar young-old median separation for the two models. Kaile v1.1 has more negative all-age and old-only Spearman point estimates, while Zihan v2.1 has higher AUC. All delta intervals include 0 except the AUC interval is bounded above by 0 because Zihan's AUC is saturated at 1.0 in this subset. This does not establish broad superiority of either model.

## Concordance Decomposition

| Subset / adjustment | Donors | Spearman | Pearson |
|---|---:|---:|---:|
| all_donors | 12 | 0.804 | 0.849 |
| young_only | 2 | NA | NA |
| old_only | 10 | 0.721 | 0.726 |
| age_group_residualized | 12 | 0.734 | 0.730 |
| age_months_residualized | 12 | 0.566 | 0.517 |

Droplet concordance remains positive within old donors and after residualizing age group, suggesting more shared donor-level structure than in FACS. However, the young-only subset has only two donors and is not interpretable.

## Technical Variable Audit

| Model | Library-size rho | Cell-count rho | Detected-gene rho | Mean cell n_counts rho | Mean cell n_genes rho |
|---|---:|---:|---:|---:|---:|
| Kaile v1.1 | -0.042 | -0.105 | 0.028 | 0.434 | 0.490 |
| Zihan v2.1 | -0.091 | -0.161 | 0.266 | 0.315 | 0.448 |

Technical correlations are lower on Droplet than on FACS, especially library-size and cell-count correlations. This supports cross-assay sensitivity, but it does not prove technical independence.

## Sex Diagnostics Among Old Donors

| Model | Old male minus female median | Wilcoxon/Mann-Whitney p-value |
|---|---:|---:|
| Kaile v1.1 | -0.064 | 0.352 |
| Zihan v2.1 | -0.081 | 0.010 |

## Droplet Review Conclusion

Both frozen FACS-derived models transport to the local Droplet limb MSC data with the expected young-old direction. Kaile v1.1 shows more negative all-age and old-only Spearman point estimates, while Zihan v2.1 reaches perfect AUC on this local subset. Paired bootstrap does not establish a broad statistically supported superiority claim. The Droplet concordance decomposition suggests shared donor-level structure beyond binary young-old separation among old donors, but this remains a sensitivity analysis rather than external validation.
