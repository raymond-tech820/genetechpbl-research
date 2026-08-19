# FACS v1.1 vs v2.1 Formal Comparison Report

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

This report compares Kaile `facs_msc_youth_score_v1_1_cleaned_limb` and Zihan `facs_msc_youth_score_v2_1` on the same cleaned TMS FACS Limb Muscle MSC input. Both models were applied as frozen scorers. No model was retrained, retuned, or modified.

## Input Data

Input file:

```text
https://drive.google.com/drive/folders/1DCOADNzh3T6XFWqdXJHKPtfXDFR_hFGl?usp=share_link
```

| Check | Value |
|---|---:|
| Cells | 815 |
| Genes | 22,966 |
| Mice | 14 |
| Diaphragm cells | 0 |
| Non-Limb_Muscle cells | 0 |
| Non-MSC cells | 0 |

## Internal Validation Evidence

This table separates validation evidence from frozen full-model same-input application.

| Model | Validation source | AUC | All-age Spearman | Old-only Spearman | Young-old median difference | Notes |
|---|---|---:|---:|---:|---:|---|
| Kaile v1.1 gene_signature | packaged donor OOF CV | 0.979 | -0.877 | -0.655 | 0.745 | OOF held-out donor scores; closest available v1.1 internal validation evidence. |
| Zihan v2.1 factorial_medium_original | nested LOMO with fold-specific feature selection | 0.979 | -0.633 | 0.655 | 0.263 | Primary predefined v2.1 validation evidence; not full-data same-input application. |

## Frozen Full-Model Same-Input Application

These results describe frozen scorer behavior and parser correctness on the same FACS donors. They should not be interpreted as independent validation.

| Model | AUC | All-age Spearman | Old-only Spearman | Young-old median difference |
|---|---:|---:|---:|---:|
| Kaile v1.1 gene_signature | 1.000 | -0.881 | -0.546 | 0.764 |
| Zihan v2.1 factorial_medium_original | 1.000 | -0.731 | 0.327 | 1.000 |

## Correctness Checks

| Check | Value |
|---|---:|
| v1.1 matched donors | 14 |
| v1.1 minimum gene overlap | 1.000 |
| v1.1 Spearman frozen application vs packaged OOF | 0.965 |
| v2.1 matched donors | 14 |
| v2.1 max absolute calibrated-score parser difference | 0.000 |
| v2.1 minimum gene coverage | 1.000 |
| v2.1 minimum weighted coverage | 1.000 |

## Paired Donor Bootstrap

Bootstrap resampled the same donors with replacement and recomputed both models' metrics on each paired resample. These intervals quantify sampling uncertainty in model differences; they are not external validation.

| Metric | Kaile observed [95% bootstrap CI] | Zihan observed [95% bootstrap CI] | Delta Kaile minus Zihan [95% bootstrap CI] |
|---|---:|---:|---:|
| AUC | 1.000 [1.000, 1.000] | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] |
| All-age Spearman | -0.881 [-0.948, -0.685] | -0.731 [-0.876, -0.268] | -0.150 [-0.527, 0.023] |
| Old-only Spearman | -0.546 [-0.907, 0.257] | 0.327 [-0.513, 0.900] | -0.873 [-1.754, 0.225] |
| Young-old median difference | 0.764 [0.531, 0.898] | 1.000 [0.866, 1.142] | -0.236 [-0.492, -0.052] |

The FACS point estimates favor Kaile v1.1 for all-age and old-only ordering, but the bootstrap intervals for the Spearman differences include 0. The young-old median difference is larger for Zihan v2.1 in same-input application.

## Concordance Decomposition

Overall model concordance can be inflated by both models separating young from old. Concordance was therefore decomposed within age groups and after residualizing age labels.

| Subset / adjustment | Donors | Spearman | Pearson |
|---|---:|---:|---:|
| all_donors | 14 | 0.684 | 0.912 |
| young_only | 6 | -0.029 | 0.148 |
| old_only | 8 | -0.286 | -0.491 |
| age_group_residualized | 14 | -0.196 | -0.291 |
| age_months_residualized | 14 | 0.222 | 0.165 |

FACS overall concordance is not preserved within young or old donors and becomes weak after age-group residualization. This means the apparent overall concordance is driven substantially by shared young-old separation.

## Technical Variable Audit

| Model | Library-size rho | Cell-count rho | Detected-gene rho | Mean cell n_counts rho | Mean cell n_genes rho |
|---|---:|---:|---:|---:|---:|
| Kaile v1.1 | -0.895 | -0.449 | 0.473 | -0.851 | 0.868 |
| Zihan v2.1 | -0.666 | -0.526 | 0.042 | -0.371 | 0.393 |

FACS technical associations are substantial for both models, especially library-size and per-cell complexity correlations. These diagnostics do not establish technical independence.

## Sex Diagnostics Among Old Donors

| Model | Old male minus female median | Wilcoxon/Mann-Whitney p-value |
|---|---:|---:|
| Kaile v1.1 | -0.241 | 0.071 |
| Zihan v2.1 | 0.140 | 0.143 |

## Feature-Level Comparison

| Module | v1.1 fold-union genes | v1.1 genes with fold frequency >= 0.6 | v2.1 genes | Same-module overlap | Same-module Jaccard | Overlap of v1.1 frequency >= 0.6 with v2.1 |
|---|---:|---:|---:|---:|---:|---:|
| young_high | 83 | 45 | 50 | 0 | 0.000 | 0 |
| old_high | 77 | 48 | 50 | 6 | 0.050 | 6 |
| any_module | 160 | 93 | 100 | 6 | 0.024 | 6 |

Gene-level overlap is very low. Pathway-level concordance is the next most informative biological comparison.

## FACS Review Conclusion

FACS supports that both frozen scorers execute correctly and separate young from old donors. Internal validation evidence should be read from Kaile packaged OOF and Zihan nested LOMO, not from full-model same-input application. Kaile v1.1 shows more consistent descriptive age ordering on FACS, especially relative to Zihan v2.1's full-model old-only direction, but paired bootstrap does not establish statistically supported superiority for Spearman-based age ordering. Neither model has demonstrated technical independence.
