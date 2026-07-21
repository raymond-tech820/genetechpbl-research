# TMS Droplet Limb Muscle MSC Youth Score v1

This folder is the GitHub-ready deliverable for the Droplet-based TMS Limb Muscle MSC Youth Score v1.

It contains the final model files, scoring code, required pseudobulk inputs for parser verification, key validation tables, figures, scripts, and the full experiment report.

## Start Here

Read:

- `Youth_score_v1_experiment_report.md`

Core final model files:

- `models/limb_msc_general_youth_score_v1_signature.csv`
- `models/limb_msc_general_youth_score_v1_calibration.json`
- `R/score_limb_msc_youth.R`

Additional exports:

- `models/limb_msc_general_youth_score_v1_signature_equal_weight_medium.csv`
- `models/limb_msc_general_youth_score_v1_large_comparator_signature.csv`

## Model Status

- Primary model: `Medium`
- Comparator model: `Large`
- Assay/source: TMS Droplet
- Tissue: `Limb_Muscle`
- Cell type: mesenchymal stem cell
- Signature size: 100 genes, 50 young-high and 50 old-high
- `bootstrap_stability = not_assessed`
- `permutation_scope = practical_training_pipeline_null`

This is a completed Droplet-based v1 modeling baseline, not yet an externally validated universal MSC aging clock.

## Quick Parser Check

From this folder:

```r
source("R/score_limb_msc_youth.R")

logcpm <- read.csv("data/pseudobulk_logcpm.csv", row.names = 1, check.names = FALSE)
scores <- score_limb_msc_youth(as.matrix(logcpm))
head(scores)
```

Expected behavior:

- 12 mouse-level pseudobulk samples are scored.
- `gene_coverage = 1`
- `weighted_coverage = 1`

Equal-weight sensitivity:

```r
scores_equal <- score_limb_msc_youth(as.matrix(logcpm), equal_weight = TRUE)
```

## Folder Layout

```text
R/          scoring function
models/     final signature, calibration JSON, comparator signatures
data/       pseudobulk logCPM/counts/metadata needed for parser checks
tables/     key Step 12-17 result tables
figures/    report figures
scripts/    scripts used to generate the workflow outputs
```

## What Is Reproducible From This Folder Alone

This folder is sufficient to:

- inspect the final v1 model and calibration;
- run the exported parser on the included pseudobulk logCPM matrix;
- verify the final Medium scores;
- review key validation, robustness, and single-cell aggregate outputs;
- inspect the scripts used for each step.

## What Requires Original Source Data

Full end-to-end reruns from raw single-cell counts require the original TMS/BPCells source data, which are not duplicated here.

In particular, full reproduction of extraction, pseudobulk construction, and single-cell scoring requires source data such as:

- original TMS droplet metadata;
- original TMS droplet BPCells matrix;
- extracted `data/limb_muscle_msc` BPCells directories.

Those files were intentionally not bundled here to keep the GitHub deliverable small and focused.

## Key Report Figures

- `figures/pseudobulk_pca_by_age.png`
- `figures/pseudobulk_pca_by_sex.png`
- `figures/volcano_plot.png`
- `figures/step15_lomo_score_by_age.png`
- `figures/step16_age_label_permutation_null.png`
- `figures/step17_single_cell_vs_pseudobulk.png`

## Main Caveats

- Young mice are only two donors, both female.
- Young-heldout folds are stress tests, not ordinary CV folds.
- The final score is sex-adjusted during feature selection, but not sex-independent.
- Cell bootstrap was not performed.
- The permutation control is practical training-pipeline permutation, not exhaustive nested outer-LOMO permutation.
- External validation on an independent dataset or different assay platform remains future work.

## Raw Single-Cell Training Data for Geneformer

The original TMS Droplet Limb Muscle MSC single-cell data have also been exported into a Geneformer-compatible AnnData format.

**Google Drive download:**
https://drive.google.com/drive/folders/1wEx3eIoCxjNavELDaVkRrRuONPeqFiqp?usp=share_link

File format:

```text
Format: AnnData (.h5ad)
Rows: individual cells
Columns: genes
Expression values: raw non-negative integer counts
Organism: Mus musculus
Tissue: Limb_Muscle
Cell type: mesenchymal stem cell
Assay: Droplet
```

The file contains:

- raw single-cell count matrix;
- mouse/donor identifiers;
- age and age group;
- sex;
- tissue and cell-type annotations;
- library metadata;
- mouse gene symbols and Ensembl identifiers, where available.

The expression matrix has not been converted to pseudobulk, logCPM, TMM-normalized expression, or z-scores. The full available gene set is retained rather than only the Youth Score signature genes.

This `.h5ad` file is intended for Geneformer or other cell-level models that require raw counts. It is different from the pseudobulk logCPM input used by the Youth Score parser.

> **Important:** The data are from *Mus musculus*. Before running Geneformer, confirm whether the selected checkpoint supports mouse genes directly or requires mouse-to-human ortholog mapping.

The training cohort has only two young mice, both female. Therefore, age and sex are partially confounded. The Youth Score feature-selection pipeline adjusts for sex, but neither the raw dataset nor the final score should be described as sex-independent.
