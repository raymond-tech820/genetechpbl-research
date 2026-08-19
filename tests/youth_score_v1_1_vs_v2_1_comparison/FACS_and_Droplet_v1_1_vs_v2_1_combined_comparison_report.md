# FACS And Droplet v1.1 vs v2.1 Combined Comparison Report

## Executive Summary

This report compares Kaile `facs_msc_youth_score_v1_1_cleaned_limb` and Zihan `facs_msc_youth_score_v2_1` on cleaned TMS FACS Limb Muscle MSC and local TMS Droplet Limb Muscle MSC. Both models were run as frozen scorers. No model was retrained, retuned, or modified.

The strongest supported conclusion is conservative: both models form directionally consistent young-old state separation on cleaned TMS FACS and local TMS Droplet donors. Kaile v1.1 shows more consistent descriptive age ordering, especially within old donors, but paired donor bootstrap does not establish statistically supported overall superiority over Zihan v2.1. Model concordance can be partly driven by shared young-old separation, especially on FACS.

FACS is same-assay application. Droplet is cross-assay sensitivity and transportability analysis. Neither analysis is independent external validation.

## Folder Organization

```text
model_comparison_v1_1_v2_1/
  facs/
  droplet/
  internal_validation_summary.csv
  paired_donor_bootstrap_summary.csv
  paired_donor_bootstrap_raw_iterations.csv
  model_concordance_decomposition.csv
  FACS_and_Droplet_v1_1_vs_v2_1_combined_comparison_report.md
```

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


## Dataset Audits

| Dataset | Cells | Genes | Mice | Diaphragm cells | Non-Limb_Muscle cells | Non-MSC cells | Notes |
|---|---:|---:|---:|---:|---:|---:|---|
| FACS | 815 | 22,966 | 14 | 0 | 0 | 0 | Diaphragm-excluded cleaned FACS limb MSC cohort |
| Droplet | 9649 | 20138 | 12 | 0 | 0 | 0 | Local subset excludes 1m and 30m donors present in packaged v1.1 Droplet sensitivity |

## Internal Validation Evidence

Validation evidence is separated from frozen full-model same-input behavior.

| Model | Validation source | AUC | All-age Spearman | Old-only Spearman | Young-old median difference | Notes |
|---|---|---:|---:|---:|---:|---|
| Kaile v1.1 gene_signature | packaged donor OOF CV | 0.979 | -0.877 | -0.655 | 0.745 | OOF held-out donor scores; closest available v1.1 internal validation evidence. |
| Zihan v2.1 factorial_medium_original | nested LOMO with fold-specific feature selection | 0.979 | -0.633 | 0.655 | 0.263 | Primary predefined v2.1 validation evidence; not full-data same-input application. |

The frozen full-data same-input FACS scores are useful for parser correctness and describing final model behavior, but they are not independent validation evidence.

## Frozen Full-Model Application And Cross-Assay Sensitivity

| Dataset | Model | AUC | All-age Spearman | Old-only Spearman | Young-old median difference |
|---|---|---:|---:|---:|---:|
| FACS same-input descriptive | Kaile v1.1 gene_signature | 1.000 | -0.881 | -0.546 | 0.764 |
| FACS same-input descriptive | Zihan v2.1 factorial_medium_original | 1.000 | -0.731 | 0.327 | 1.000 |
| Droplet cross-assay sensitivity | Kaile v1.1 gene_signature | 0.950 | -0.873 | -0.817 | 0.126 |
| Droplet cross-assay sensitivity | Zihan v2.1 factorial_medium_original | 1.000 | -0.830 | -0.701 | 0.125 |

## Paired Donor Bootstrap Difference Analysis

Each bootstrap iteration resampled the same donors with replacement and recomputed both models' metrics. Delta is Kaile minus Zihan.

### FACS

| Metric | Kaile observed [95% bootstrap CI] | Zihan observed [95% bootstrap CI] | Delta Kaile minus Zihan [95% bootstrap CI] |
|---|---:|---:|---:|
| AUC | 1.000 [1.000, 1.000] | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] |
| All-age Spearman | -0.881 [-0.948, -0.685] | -0.731 [-0.876, -0.268] | -0.150 [-0.527, 0.023] |
| Old-only Spearman | -0.546 [-0.907, 0.257] | 0.327 [-0.513, 0.900] | -0.873 [-1.754, 0.225] |
| Young-old median difference | 0.764 [0.531, 0.898] | 1.000 [0.866, 1.142] | -0.236 [-0.492, -0.052] |

### Droplet

| Metric | Kaile observed [95% bootstrap CI] | Zihan observed [95% bootstrap CI] | Delta Kaile minus Zihan [95% bootstrap CI] |
|---|---:|---:|---:|
| AUC | 0.950 [0.727, 1.000] | 1.000 [1.000, 1.000] | -0.050 [-0.273, 0.000] |
| All-age Spearman | -0.873 [-0.964, -0.624] | -0.830 [-0.954, -0.384] | -0.044 [-0.412, 0.208] |
| Old-only Spearman | -0.817 [-0.953, -0.427] | -0.701 [-0.895, -0.037] | -0.117 [-0.691, 0.342] |
| Young-old median difference | 0.126 [0.071, 0.176] | 0.125 [0.058, 0.174] | 0.002 [-0.033, 0.060] |

Interpretation: Kaile v1.1 has more negative age-ordering point estimates in several comparisons, but the bootstrap intervals for Spearman differences include 0. Therefore the current data do not prove Kaile v1.1 is statistically superior overall. Zihan v2.1 has a larger same-input FACS young-old median difference and saturated Droplet AUC, but these also need cautious interpretation.

## Concordance Decomposition

Overall model concordance can be inflated when both models separate young from old. Concordance was decomposed within age strata and after residualizing age labels.

### FACS Concordance

| Subset / adjustment | Donors | Spearman | Pearson |
|---|---:|---:|---:|
| all_donors | 14 | 0.684 | 0.912 |
| young_only | 6 | -0.029 | 0.148 |
| old_only | 8 | -0.286 | -0.491 |
| age_group_residualized | 14 | -0.196 | -0.291 |
| age_months_residualized | 14 | 0.222 | 0.165 |

### Droplet Concordance

| Subset / adjustment | Donors | Spearman | Pearson |
|---|---:|---:|---:|
| all_donors | 12 | 0.804 | 0.849 |
| young_only | 2 | NA | NA |
| old_only | 10 | 0.721 | 0.726 |
| age_group_residualized | 12 | 0.734 | 0.730 |
| age_months_residualized | 12 | 0.566 | 0.517 |

FACS overall concordance is largely not preserved within age groups or after age-group residualization, so it should not be interpreted as full biological agreement. Droplet concordance remains positive among old donors and after age-group residualization, suggesting stronger shared donor-level structure in the cross-assay setting, though this is still a small-donor sensitivity result.

## Technical Variable Audit

| Dataset | Model | Library-size rho | Cell-count rho | Detected-gene rho | Mean cell n_counts rho | Mean cell n_genes rho |
|---|---|---:|---:|---:|---:|---:|
| FACS | Kaile v1.1 | -0.895 | -0.449 | 0.473 | -0.851 | 0.868 |
| FACS | Zihan v2.1 | -0.666 | -0.526 | 0.042 | -0.371 | 0.393 |
| Droplet | Kaile v1.1 | -0.042 | -0.105 | 0.028 | 0.434 | 0.490 |
| Droplet | Zihan v2.1 | -0.091 | -0.161 | 0.266 | 0.315 | 0.448 |

Technical correlations are much stronger on FACS than on Droplet. This does not prove technical independence on Droplet, but it suggests that the strong FACS library-size axis is not identically reproduced in the Droplet sensitivity setting.

## Sex Diagnostics Among Old Donors

| Dataset | Model | Old male minus female median | Wilcoxon/Mann-Whitney p-value |
|---|---|---:|---:|
| FACS | Kaile v1.1 | -0.241 | 0.071 |
| FACS | Zihan v2.1 | 0.140 | 0.143 |
| Droplet | Kaile v1.1 | -0.064 | 0.352 |
| Droplet | Zihan v2.1 | -0.081 | 0.010 |

## Feature-Level Comparison On FACS

| Module | v1.1 fold-union genes | v1.1 genes with fold frequency >= 0.6 | v2.1 genes | Same-module overlap | Same-module Jaccard | Overlap of v1.1 frequency >= 0.6 with v2.1 |
|---|---:|---:|---:|---:|---:|---:|
| young_high | 83 | 45 | 50 | 0 | 0.000 | 0 |
| old_high | 77 | 48 | 50 | 6 | 0.050 | 6 |
| any_module | 160 | 93 | 100 | 6 | 0.024 | 6 |

Gene-level overlap is very low, especially for young-high genes. The two models can agree at the donor/ranking level while using largely distinct gene sets. This makes pathway-level concordance the most informative next biological comparison.

## Overall Comparison Conclusion

The current data support this conservative conclusion:

Both models can be executed correctly and both form directionally consistent young-old state separation on cleaned TMS FACS and local TMS Droplet donors. Droplet results show cross-assay directional transportability for both models.

The current data do not support a strong claim that Kaile v1.1 has been proven superior to Zihan v2.1. Kaile v1.1 shows more consistent descriptive age ordering, especially within old donors, but paired donor bootstrap intervals do not establish statistically supported superiority for Spearman-based age-ordering metrics.

Zihan v2.1 remains a valid frozen cohort-state separator with strong young-old separation, but its FACS old-only behavior argues against interpreting it as a continuous aging trajectory score.

Neither model has demonstrated technical independence. Neither model should be described as externally validated from these analyses.

The low gene-level overlap and moderate-to-high donor-level concordance make pathway-level concordance the next most useful biological analysis.

## Report Files

```text
facs/FACS_v1_1_vs_v2_1_formal_comparison_report.md
droplet/Droplet_v1_1_vs_v2_1_formal_comparison_report.md
FACS_and_Droplet_v1_1_vs_v2_1_combined_comparison_report.md
internal_validation_summary.csv
paired_donor_bootstrap_summary.csv
model_concordance_decomposition.csv
```
