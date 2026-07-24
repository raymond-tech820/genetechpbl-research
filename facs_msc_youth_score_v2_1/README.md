# FACS MSC Youth Score v2.1

Frozen diaphragm-excluded TMS FACS Limb Muscle MSC Youth Score model.

## Model Status

This folder contains the minimal files needed to apply the **FACS v2.1** frozen model. The primary model is:

`factorial_medium_original`

Interpretation boundary: this model is a reproducible **young-old cohort-state separator**. It is **not** a validated continuous aging clock, and technical independence has not been demonstrated.

## Training Data

FACS v2.1 was trained from the diaphragm-excluded TMS FACS Limb Muscle MSC raw-count training set. The cleaned training data contained:

- 815 cells
- 22,966 genes
- 14 mice
- 264 young 3m cells
- 551 old 18m/24m cells
- 0 diaphragm cells

All remaining `subtissue` values are forelimb/hindlimb formatting variants, standardized during analysis as `Forelimb_Hindlimb`.

## Folder Contents

- `models/`: frozen signatures, calibration parameters, and individual signature CSVs.
- `parser/score_facs_v2_1_youth_model.R`: R parser for gene x sample pseudobulk raw counts.
- `parser/score_facs_v2_1_single_cell_h5ad.py`: Python wrapper for single-cell raw-count h5ad input; internally pseudobulks by donor/sample.
- `data/`: small pseudobulk training example and frozen training scores for parser checks. This is not raw single-cell data.
- `tables/`: key result summaries used by the report.
- `figures/`: selected figures used by the report.
- `reports/`: step-level reports and model card.
- `FACS_Youth_Score_v2_1_experiment_report.md`: final experiment report.


## Original Training Data

Raw single-cell training data are not included in this GitHub folder. Add the download location here before sharing with teammates who need to rerun preprocessing or retraining.

Training data URL: `https://drive.google.com/drive/folders/1DCOADNzh3T6XFWqdXJHKPtfXDFR_hFGl?usp=share_link`

## Apply To Pseudobulk Counts In R

Input should be a raw-count matrix with genes as rows and samples/donors as columns. Counts should not be CPM-normalized, log-transformed, z-scored, or TMM-normalized.

```r
source("parser/score_facs_v2_1_youth_model.R")

counts_df <- read.csv("data/facs_v2_1_limb_msc_pseudobulk_counts.csv", check.names = FALSE)
counts <- as.matrix(counts_df[, -1])
rownames(counts) <- counts_df$gene

scores <- score_facs_v2_1_youth_model(
  counts = counts,
  signature_csv = "models/facs_v2_1_full_data_frozen_signatures_all_models.csv",
  calibration_csv = "models/facs_v2_1_full_data_frozen_calibration.csv",
  model = "factorial_medium_original"
)
print(scores)
```

## Apply To Single-cell h5ad In Python

Input h5ad must contain raw non-negative integer-like counts in `X`, `raw.X`, or a named layer. The matrix orientation must be cells x genes. The sample column should identify the donor or mouse to pseudobulk.

```bash
python parser/score_facs_v2_1_single_cell_h5ad.py   --input-h5ad path/to/raw_counts.h5ad   --sample-column mouse.id   --signature-csv models/facs_v2_1_full_data_frozen_signatures_all_models.csv   --calibration-csv models/facs_v2_1_full_data_frozen_calibration.csv   --model factorial_medium_original   --output scores/facs_v2_1_scores.csv
```

## Score Formula

For each input sample, the parser computes deployable expression:

```text
E_g = log2(CPM_g + 1)
Z_g = (E_g - mu_g) / s_g
```

where `mu_g` and `s_g` are frozen training-set mean and standard deviation stored in the signature file.

Module scores are weighted means:

```text
M_young = sum(w_g Z_g for young-high genes) / sum(w_g)
M_old   = sum(w_g Z_g for old-high genes) / sum(w_g)
S_raw   = M_young - M_old
S_cal   = (S_raw - C_old) / (C_young - C_old)
```

Higher calibrated scores indicate a more young-like transcriptional state relative to the FACS v2.1 training medians.

Because calibration uses training medians, full-data training-set medians are fixed by construction. Validation evidence should be read from nested LOMO, bootstrap, and permutation results, not apparent full-data separation.

## Key Internal Results

Primary nested LOMO model, `factorial_medium_original`:

- AUC young vs old: 0.979
- all-age Spearman: -0.633
- old-only Spearman: 0.655
- library-size Spearman: -0.521
- zero-signature folds: 0

Formal 999 nested permutation:

- primary statistic `abs_all_age_rho`: observed 0.633, empirical p = 0.153
- supporting `abs(AUC - 0.5)`: observed 0.479, empirical p = 0.039
- recorded young-minus-old median: empirical p = 0.002

Conclusion: young-old separation has moderate support, but the pre-locked primary permutation test does not support a robust continuous all-age trajectory beyond the practical randomized-label null.

## Recommended Reporting Language

Use:

> The diaphragm-excluded FACS v2.1 model is a frozen FACS Limb Muscle MSC cohort-state score with reproducible young-old separation.

Do not use:

> Validated aging clock.

> Technical-confounding-free model.

> Externally validated model.

## Dependencies

R parser:

- Base R only for scoring from pseudobulk counts.

Python h5ad wrapper:

- `anndata`
- `numpy`
- `pandas`
- `scipy`

## Model Versions

Available model names:

- `factorial_medium_original` - primary
- `factorial_large_original` - size comparator
- `factorial_stability_selected` - stability comparator
- `factorial_medium_equal_weight` - weight sensitivity
- `age_only_de` - age-only baseline
