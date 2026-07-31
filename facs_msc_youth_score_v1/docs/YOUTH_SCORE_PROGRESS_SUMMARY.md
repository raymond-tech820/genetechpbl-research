# Youth Score Project Progress Summary

**Updated: 2026-07-17**  
**Scope: TMS MSC age scoring, cross-modality sensitivity analysis, and GSE176206 external validation**

## 1. Executive summary

This phase used two tissue-specific MSC populations from the local Tabula Muris Senis (TMS) FACS data to develop Youth Scores: SCAT adipose MSC as the primary analysis and Limb Muscle MSC as the secondary analysis. We compared three candidate models—Gene signature, Elastic Net, and a Gene-token Transformer—and used a Technical-only logistic regression as a diagnostic control for technical confounding.

All internal performance estimates are based on out-of-sample predictions from **mouse-grouped, age-stratified five-fold cross-validation**, not training-set fit. Elastic Net was selected for SCAT (donor-level ROC-AUC 0.946), while Gene signature was selected for Limb Muscle (ROC-AUC 1.000). The Transformer learned an age signal in both populations but did not outperform the simpler models.

When the Limb FACS model was applied to TMS Droplet Limb MSC, donor-level Youth Score was strongly negatively correlated with age (Spearman \(\rho=-0.856\)). This supports directional transfer across assay technologies, although age, sex, donor count, and modality remain partially confounded.

In the training-naive GSE176206 data, mean Youth Score was higher in Young control than Aged control for both Adipo and MSC. The MSC age contrast was +0.030 with a 95% donor-bootstrap confidence interval of +0.0135 to +0.0453. Adipo had only one pooled library per condition in the processed metadata, preventing a donor-level confidence interval. The external SOKM analysis did not show movement of aged samples toward the current TMS-defined youth axis.

## 2. Data and analysis design

### 2.1 Formal training data

SCAT and Limb Muscle are not equal halves of one MSC dataset. They are two independently selected, tissue-specific MSC populations from the FACS data.

| Analysis population | Assay | Cells | Young mice | Old mice | Age definition | Role |
|---|---:|---:|---:|---:|---|---|
| SCAT adipose MSC | FACS | 1,570 | 7 | 8 | Young=3m; Old=18m/24m | Formal primary analysis |
| Limb Muscle MSC | FACS | 935 | 6 | 8 | Young=3m; Old=18m/24m | Formal secondary analysis |

The training label is Young/Old, and Youth Score is defined as fold-calibrated `P(young)`. Under external domain shift, the score is used primarily for relative comparisons within the same external dataset and should not be interpreted directly as a biological-age probability.

### 2.2 Four analysis strategies

| Method | Category | Role in this project |
|---|---|---|
| Gene signature | Bioinformatic score + logistic regression | Uses training mice only to identify youth- and age-associated genes and construct a one-dimensional score |
| Elastic Net | Regularized logistic regression | Learns sparse linear weights from gene expression |
| Gene-token Transformer | Deep-learning Transformer | Learns nonlinear combinations of gene identity, expression bins, and expression rank, with auxiliary age prediction |
| Technical only | Classical logistic regression | Uses total counts, detected genes, and sex to diagnose technical confounding |

The first three are formal candidate models. Technical only is a diagnostic control and is **not eligible to become the official model**.

### 2.3 Exact meaning of five-fold cross-validation

We used **mouse-grouped, age-stratified five-fold outer cross-validation with an inner validation split from the outer training mice**:

1. `mouse.id` is the minimum splitting unit, so cells from one mouse cannot appear in different partitions.
2. The 15 SCAT mice or 14 Limb mice are divided into five outer folds, each taking one turn as the independent test fold.
3. From the remaining mice, at least one Young and one Old mouse are reserved for early stopping, probability calibration, and threshold selection.
4. HVG selection, normalization, gene signatures, expression bins, and probability calibration use only the relevant training/validation mice.
5. At the end of five folds, every mouse has exactly one fully out-of-fold prediction; the final AUC pools these donor-level out-of-fold predictions.

This is neither a one-time 50/50 train-validation split nor a random cell-level split.

## 3. Internal validation in TMS FACS

| Method | SCAT AUC | Limb Muscle AUC |
|---|---:|---:|
| Gene signature | 0.786 | **1.000** |
| Elastic Net | **0.946** | 0.979 |
| Transformer | 0.857 | 0.875 |
| Technical only | 0.696 | 0.688 |

![Internal donor-level ROC-AUC comparison](../outputs/youth_score/progress_summary/figures/figure_1_internal_auc.png)

*Figure 1. Points are pooled out-of-fold donor-level ROC-AUC values; horizontal bars are 95% intervals from 2,000 donor-bootstrap iterations; the dashed line marks chance at 0.5; black-edged points denote the selected models. Vector version: [PDF](../outputs/youth_score/progress_summary/figures/figure_1_internal_auc.pdf).*

### 3.1 Interpretation

- The official SCAT model is Elastic Net: ROC-AUC 0.946, 95% interval 0.821–1.000, balanced accuracy 0.804, and Brier score 0.088.
- The official Limb Muscle model is Gene signature: ROC-AUC 1.000, balanced accuracy 1.000, and Brier score 0.047.
- All three gene-expression candidate models have higher AUC point estimates than Technical only in both populations, indicating that the formal models use information beyond sequencing depth, detected gene count, and sex.
- With only 15 and 14 donors, several AUC intervals are broad and overlap. A higher point estimate should therefore not be described as proof that every candidate is statistically significantly better than Technical only.
- The Transformer showed no performance advantage. Many cells do not replace independent mice, and the number of independent donors may be too small to favor a complex neural network.
- Limb AUC=1.000 means perfect ranking of the current 14 held-out mice, not guaranteed 100% accuracy in future data.

## 4. TMS Droplet cross-modality sensitivity analysis

The official Limb FACS Gene-signature model was frozen and applied directly to separately processed TMS Droplet Limb Muscle MSC. Droplet data were not used for formal training, feature selection, or calibration.

| Age | Mean Youth Score | Mice |
|---|---:|---:|
| 1m | 0.448 | 2 |
| 3m | 0.273 | 2 |
| 18m | 0.215 | 4 |
| 21m | 0.162 | 2 |
| 24m | 0.112 | 4 |
| 30m | 0.142 | 2 |

The donor-level Spearman correlation is:

```text
rho = -0.856
```

![Droplet age trajectory](../outputs/youth_score/progress_summary/figures/figure_2_droplet_age_trajectory.png)

*Figure 2. Small points are individual mice; open points and the connecting line are unweighted mouse means for each age. Vector version: [PDF](../outputs/youth_score/progress_summary/figures/figure_2_droplet_age_trajectory.pdf).*

### 4.1 Interpretation

Youth Score generally decreases with age, as expected, although 30m is slightly higher than 24m and the trajectory is not strictly monotonic. This supports preservation of the age direction after transfer from FACS to Droplet. It should be described as a **within-TMS cross-modality sensitivity analysis**, not a fully independent cross-dataset validation. Sex composition and mouse counts differ by age, and only 2–4 mice are available per age, so the correlation alone cannot establish a stable biological-age effect.

## 5. External age-direction validation in GSE176206

[GSE176206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176206) is an independent mouse partial-reprogramming single-cell RNA-seq study containing Young/Aged and Control/SOKM adipogenic cells and limb muscle MSCs. External data were not used for training, feature selection, or calibration.

The model-to-dataset mapping is:

| TMS model | GSE176206 data | Match quality |
|---|---|---|
| SCAT adipose MSC / Elastic Net | Adipo SOKM | Adipose tissue and lineage related, but not an identical MSC identity; an approximate cross-cell-state validation |
| Limb Muscle MSC / Gene signature | MSC SOKM | Limb-muscle-derived MSC in both datasets; a more direct match |

### 5.1 Young control versus Aged control

| External data | Young control | Aged control | Young − Aged |
|---|---:|---:|---:|
| Adipo | 0.497 | 0.447 | +0.050 |
| MSC | 0.493 | 0.463 | +0.030 |

![External age-direction validation](../outputs/youth_score/progress_summary/figures/figure_3_external_age_validation.png)

*Figure 3. The left panel shows condition scores and available donor points; the right panel shows the Young control − Aged control effect. The MSC bar is a 95% donor-bootstrap interval. The open Adipo diamond indicates that the processed metadata contain only one pooled library per condition, so no donor-level interval is available. Vector version: [PDF](../outputs/youth_score/progress_summary/figures/figure_3_external_age_validation.pdf).*

### 5.2 Correct interpretation of the confidence interval

The MSC contrast is:

```text
Young control - Aged control = +0.02985
95% donor-bootstrap CI       = [+0.01350, +0.04532]
```

`+0.0135 to +0.0453` is the uncertainty interval for the **MSC contrast estimate of +0.030 itself**. The point estimate lies inside this interval. It is incorrect to say that the Adipo or MSC differences are partly or entirely greater than this interval, and the Adipo difference should not be compared directly with the MSC confidence interval. Because the entire MSC interval is above zero, the data support a higher score in MSC Young control than Aged control.

Adipo has one pooled library per condition in the processed metadata. Thousands of cells do not replace independent mice, so +0.050 is a directional observation without a reliable donor-level confidence interval.

### 5.3 Supported conclusion

The age direction is correct in both external datasets, providing preliminary evidence of cross-study transfer. This validates the **relative age direction**, not external classification accuracy: no external accuracy or ROC-AUC was calculated from these means, and platform and cell-state shifts make a fixed 0.5 threshold inappropriate.

MSC is the stronger and better cell-type-matched external result. Adipo is a weaker cross-cell-state check.

## 6. GSE176206 SOKM analysis

The SOKM effect is defined as:

```text
Aged SOKM - Aged control
```

Only a positive value would indicate movement toward the current TMS-defined youth axis.

| Data | Aged control | Aged SOKM | SOKM − Control |
|---|---:|---:|---:|
| Adipo | 0.447 | 0.436 | −0.011 |
| MSC | 0.463 | 0.344 | −0.118 |

![External SOKM contrast](../outputs/youth_score/progress_summary/figures/figure_4_external_sokm_contrast.png)

*Figure 4. The left panel shows aged-control and aged-SOKM scores; the right panel shows Aged SOKM − Aged control. The MSC bar is a 95% donor-bootstrap interval; Adipo has no donor-resolved interval. Vector version: [PDF](../outputs/youth_score/progress_summary/figures/figure_4_external_sokm_contrast.pdf).*

The MSC treatment contrast is −0.118 with a 95% donor-bootstrap interval of −0.161 to −0.055, entirely below zero. The Adipo contrast is −0.011 but lacks a usable independent-donor interval.

The strict conclusion supported by the current analysis is:

> In GSE176206, SOKM did not move aged samples toward the youth axis learned by this project from natural aging in TMS; the observed MSC shift was in the opposite direction.

This does not prove that SOKM accelerates aging or is ineffective on other endpoints. Potential explanations include differences between natural-aging and reprogramming-response axes, acute SOKM stress, temporary suppression of cell identity, assay and culture domain shifts, imperfect Adipo-to-SCAT-MSC identity matching, and donor/cell-count imbalance. The original study reported youthful movement using its own age-expression programs; the present result asks whether treatment aligns with an **independently trained TMS Youth Score**, which is a different statistical question. [Roux et al., 2022](https://doi.org/10.1016/j.cels.2022.05.002)

## 7. Graded conclusions

| Scientific question | Current result | Evidence strength |
|---|---|---|
| Can Young and Old MSCs be distinguished in TMS FACS? | Supported; SCAT AUC 0.946 and Limb AUC 1.000 | Strong within this study, but donor-limited |
| Does the Transformer outperform simpler models? | Not supported | Clear |
| Does the Limb FACS model retain an age direction in Droplet? | Supported, \(\rho=-0.856\) | Moderate; confounded |
| Is external MSC Young control higher than Aged control? | Supported, +0.030 with a 95% interval above zero | Preliminary; three animals per condition |
| Is the age direction correct in external Adipo? | Directionally yes, +0.050 | Weak; one pooled library per condition and approximate cell identity |
| Does SOKM increase Youth Score in aged samples? | Not supported; opposite direction in MSC | Clearer for MSC, uncertain for Adipo |

## 8. Claims that are and are not supported

### Supported statements

- We developed and validated two tissue-specific TMS MSC Youth Scores.
- Simpler models outperformed the Transformer at the current donor scale.
- The Limb model retained a general age-decreasing trajectory in TMS Droplet.
- External MSC Young control scored higher than Aged control.
- The current model did not detect an SOKM-driven increase in Youth Score in GSE176206.

### Unsupported statements

- Youth Score cannot be interpreted as clinical age, treatment safety, causal rejuvenation, or functional recovery.
- Thousands of cells cannot be treated as thousands of biological replicates.
- Limb AUC=1.000 does not guarantee 100% accuracy in future data.
- A negative SOKM contrast does not prove that SOKM is harmful or accelerates aging.
- The Adipo result is not a fully cell-type-matched external validation of SCAT MSC.

## 9. Files and reproducibility

All scientific figures are available as 300 dpi PNG and vector PDF files. They are generated automatically from existing result tables by the following English-language script:

```powershell
.\.venv\Scripts\python.exe scripts\plot_youth_score_progress.py
```

- Plotting script: [`scripts/plot_youth_score_progress.py`](../scripts/plot_youth_score_progress.py)
- Figure directory: [`outputs/youth_score/progress_summary/figures`](../outputs/youth_score/progress_summary/figures)
- SCAT metrics: [`outputs/youth_score/scat_facs/model_metrics.csv`](../outputs/youth_score/scat_facs/model_metrics.csv)
- Limb metrics: [`outputs/youth_score/limb_facs/model_metrics.csv`](../outputs/youth_score/limb_facs/model_metrics.csv)
- Droplet donor results: [`outputs/youth_score/limb_droplet_sensitivity/donor_scores.csv`](../outputs/youth_score/limb_droplet_sensitivity/donor_scores.csv)
- External summary: [`outputs/youth_score/external/GSE176206/combined_summary.json`](../outputs/youth_score/external/GSE176206/combined_summary.json)

