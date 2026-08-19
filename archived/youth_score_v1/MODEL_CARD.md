# Model Card: FACS MSC Youth Score v1

## Purpose

This release packages two research-only transcriptional Youth Scores trained on TMS FACS mesenchymal stem-cell (MSC) cohorts. A higher score indicates a more young-like transcriptomic profile relative to the corresponding 3-month TMS reference.

## Reference data and labels

| Model | Tissue and cell type | Young label | Old label | Cells | Young / old donors |
| --- | --- | --- | --- | ---: | ---: |
| `scat_facs` | SCAT, mesenchymal stem cell of adipose | 3 months | 18 and 24 months | 1,570 | 7 / 8 |
| `limb_facs` | Limb Muscle, mesenchymal stem cell | 3 months | 18 and 24 months | 935 | 6 / 8 |

The source matrix is TMS FACS raw count data. Official TMS QC was retained; zero-count cells and cells with fewer than 500 matrix-detected genes were removed.

## Modeling and selection

Each cohort was evaluated with five donor-grouped outer folds. All cells from a mouse stayed in the same fold. Each training split used only its training donors for HVG selection, tokenization, calibration, and signature construction.

Four methods were evaluated: gene signature, Elastic Net logistic regression, a compact gene-token Transformer with an auxiliary age head, and a technical-only diagnostic model based on total counts, detected genes, and sex. The official model was selected by donor-level out-of-fold ROC-AUC; AUC differences below 0.02 were resolved by lower Brier score, then by simplicity.

| Method | SCAT ROC-AUC | Limb ROC-AUC |
| --- | ---: | ---: |
| Gene signature | 0.786 | 1.000 |
| Elastic Net | 0.946 | 0.979 |
| Gene Transformer | 0.857 | 0.875 |
| Technical-only diagnostic | 0.696 | 0.688 |

Therefore, the official SCAT model is Elastic Net and the official Limb model is the gene-signature score. Transformer checkpoints are included for reproducibility and auxiliary age prediction, but the Transformer was not the selected Youth Score model.

## Input and output

Input must be a cell-by-gene raw-count matrix. The CLI reads `.raw.X` when present in an h5ad file, otherwise `.X`. Internally it applies cell-level library-size normalization to 10,000 and `log1p`; input must not be pre-normalized or log-transformed. Gene matching is symbol based. A minimum 70% overlap with every fold vocabulary is required.

The score is the mean calibrated `P(young)` across the five fold-specific models. `predicted_age_months` is an auxiliary Transformer output and must not be interpreted as an accurate chronological-age estimator.

## External checks

The FACS Limb model showed an age-ordered trend when applied to TMS Droplet Limb MSC donors (rho = -0.856, 16 mice). The GSE176206 Adipo and MSC young-control groups both scored above their matched aged-control groups. The MSC young-minus-aged-control difference was +0.0298 with donor bootstrap interval [+0.0135, +0.0453].

For aged SOKM versus aged control, scores changed by -0.0109 (Adipo; one pooled library per condition) and -0.1185 (MSC; 95% donor interval [-0.1612, -0.0551]). These data do not support a young-direction shift under this score.

## Limitations and non-claims

- This is not a clinical age predictor, safety score, treatment recommendation, or causal rejuvenation measure.
- Donor numbers are limited. Age, sex, assay, and technical covariates are partly confounded in TMS.
- Cells are not biological replicates. Aggregate scores by donor or true replicate before inferential comparisons.
- Cross-dataset absolute scores are not directly interchangeable. Prefer within-dataset contrasts and check gene overlap, cell identity, and reference-expression compatibility.
- The reported AUCs are internal donor-level out-of-fold estimates, not proof of general clinical performance.
