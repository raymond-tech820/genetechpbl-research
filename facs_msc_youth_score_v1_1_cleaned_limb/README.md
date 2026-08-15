# FACS Limb Muscle MSC Youth Score v1.1 (cleaned cohort)

This release updates the FACS Limb Muscle MSC component of the prior FACS Youth Score deliverable after a metadata audit identified 120 old cells from `Muscle Diaphragm` that were coarsely labelled as `Limb_Muscle`.

## What changed

- Excluded all cells with `subtissue = Muscle Diaphragm` before feature selection, calibration, and training.
- Retrained all five donor-grouped outer folds and all three Transformer seeds from the corrected raw-count bundle.
- Re-ran the primary analysis, all-gene sensitivity analysis, FACS-to-Droplet sensitivity check, and GSE176206 MSC external validation.
- Preserved the prior release rather than overwriting it; this directory is the corrected replacement for its Limb component.

## Corrected reference cohort

| Item | Value |
| --- | ---: |
| Assay | TMS FACS |
| Tissue / cell type | Limb_Muscle / mesenchymal stem cell |
| Young label | 3 months |
| Old label | 18 and 24 months |
| Cells | 815 (264 young, 551 old) |
| Biological donors | 14 (6 young, 8 old) |
| Excluded cells | 120 old diaphragm cells |

## Primary donor-level cross-validation

| Method | ROC-AUC | Brier score |
| --- | ---: | ---: |
| Gene-signature | **0.979** | **0.055** |
| Elastic Net | 0.979 | 0.072 |
| Gene-token Transformer | 0.854 | 0.185 |
| Technical-only diagnostic | 0.583 | 0.272 |

The official corrected scorer is the gene-signature model: it tied for highest donor-level ROC-AUC and had the lower Brier score. The score is a research-only transcriptomic similarity measure to the TMS 3-month reference, not a clinical age estimate or safety claim.

## Independent checks

- **FACS-to-Droplet Limb MSC:** donor-level age correlation `Spearman rho = -0.838`; this supports the expected age direction across platforms but remains a sensitivity analysis because age, sex, donor count, and platform are partly confounded.
- **GSE176206 MSC:** Young control scored above Aged control by `+0.036` (95% donor-bootstrap interval `[+0.019, +0.052]`).
- **GSE176206 Aged SOKM:** Aged SOKM scored below Aged control by `-0.149` (95% donor-bootstrap interval `[-0.197, -0.077]`). Under this score, the data do not support a rejuvenation-like SOKM shift. This is not evidence that SOKM is broadly harmful or accelerates aging.

## Contents

- `models/limb_facs/`: fold-specific models, tokenizers, calibration, out-of-fold predictions, metrics, and all-gene sensitivity results.
- `results/`: corrected FACS-to-Droplet and GSE176206 MSC validation outputs.
- `reports/`: generated training and external-validation reports.
- `research_pipeline/`: reproducible Python package, configurations, dependency lock file, and relevant tests.

Raw TMS and GSE176206 expression matrices are not redistributed here. They must be obtained from their original public sources and placed according to the paths documented in `research_pipeline/configs/`.
