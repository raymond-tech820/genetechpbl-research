# TMS FACS Limb Muscle MSC Youth Score v2

This folder is the GitHub-ready deliverable for the frozen **FACS-derived Limb Muscle MSC Youth Score v2**.

The primary model is `factorial_medium_original`: 50 young-high and 50 old-high genes trained on mouse-level pseudobulk from the TMS FACS Limb_Muscle MSC subset. Comparator models are included for audit and sensitivity, but the primary model role is frozen.

## What Is Included

- `models/`: frozen signatures, calibration parameters, and RDS bundle.
- `parser/`: deployable R scoring function.
- `data/`: compact metadata and frozen full-data training scores; large training count matrices are intentionally omitted from the GitHub package.
- `tables/`: compact validation and summary tables.
- `figures/`: selected QC, validation, and cross-assay figures used by the report.
- `reports/`: intermediate reports and model card.
- `scripts/`: scripts used to reproduce the v2 analysis chain.
- `FACS_Youth_Score_v2_experiment_report.md`: final experiment report.

## Raw Data Scope

The GitHub package does not include raw-count matrices. The filtered Limb_Muscle mesenchymal stem cell subset used for v2 training is stored separately as the Geneformer/raw-data bundle.

- Genes: 22966
- Cells: 935
- Mice: 14
- Young cells: 264
- Old cells: 671
- Matrix/metadata cell order match: TRUE

Geneformer-oriented raw-data files are intentionally not included in this GitHub deliverable. They are stored separately for Google Drive upload.

Google Drive URL: `https://drive.google.com/drive/folders/1q2xSg7y7pWaGOFbT6cpl7_wFkNQeL2q7?usp=share_link`

Expected Geneformer raw-data files in the Google Drive bundle:

- `tms_facs_limb_msc_raw_counts.h5ad`: AnnData, cells x genes.
- `tms_facs_limb_msc_raw_counts_genes_by_cells.mtx`: Matrix Market, genes x cells.
- `tms_facs_limb_msc_cell_metadata.csv`
- `tms_facs_limb_msc_gene_metadata.csv`
- `tms_facs_limb_msc_genes.txt`
- `tms_facs_limb_msc_cells.txt`

## Requirements

For applying the frozen model, only base R is required by the parser. The input must be a genes x samples raw-count pseudobulk matrix with gene symbols as row names. Single cells should first be aggregated to donor/sample-level pseudobulk before model-level interpretation.

The training/audit scripts are included for transparency. They expect the original `Youth_score` project layout plus the FACS/Droplet data inputs or the external Geneformer/raw-data bundle; they are not required for routine scoring with the frozen parser.

## Minimal Scoring Example

```r
source("parser/score_facs_v2_youth_model.R")
# Provide a genes x samples pseudobulk raw-count matrix with gene symbols as rownames.
# For example, build this from the Google Drive Geneformer raw-data bundle or an external dataset.
counts <- your_pseudobulk_count_matrix
score <- score_facs_v2_youth_model(
  counts,
  signature_csv = "models/facs_v2_full_data_frozen_signatures_all_models.csv",
  calibration_csv = "models/facs_v2_full_data_frozen_calibration.csv",
  model = "factorial_medium_original"
)
head(score)
# Check gene_coverage and weighted_coverage before interpreting scores.
```

## Interpretation Boundary

This is an internally finalized and frozen model with within-TMS cross-assay sensitivity evidence. It is **not** an externally validated universal MSC aging clock. Technical and cohort-structure limitations are documented in the report and model card.
