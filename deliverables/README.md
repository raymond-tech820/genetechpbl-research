Limb Muscle MSC Youth Score v1

The complete Droplet-based Limb Muscle MSC Youth Score model is available in:

limb_msc_youth_score_v1/

This folder contains:

* the finalized 100-gene Medium Youth Score signature;
* the Large high-stability comparator;
* the equal-weight Medium sensitivity version;
* calibration parameters;
* the R scoring parser;
* pseudobulk input data;
* validation and robustness results;
* figures, tables, and reproducible scripts for Steps 1–17.

The primary parser takes a gene-by-sample pseudobulk logCPM matrix and returns:

* raw Youth Score;
* calibrated Youth Score;
* clipped score;
* gene coverage;
* weighted gene coverage.

See the experiment report inside the folder for the complete methodology, validation results, limitations, and usage instructions.

Raw Single-Cell Training Data for Geneformer

The original TMS Droplet Limb Muscle MSC single-cell data have also been exported into a Geneformer-compatible AnnData format.

Google Drive download:
[
Download tms_droplet_limb_msc_raw_counts.h5ad](https://drive.google.com/drive/folders/1wEx3eIoCxjNavELDaVkRrRuONPeqFiqp?usp=share_link)

File format:

Format: AnnData (.h5ad)
Rows: individual cells
Columns: genes
Expression values: raw non-negative integer counts
Organism: Mus musculus
Tissue: Limb_Muscle
Cell type: mesenchymal stem cell
Assay: Droplet

The file contains:

* raw single-cell count matrix;
* mouse/donor identifiers;
* age and age group;
* sex;
* tissue and cell-type annotations;
* library metadata;
* mouse gene symbols and Ensembl identifiers, where available.

The expression matrix has not been converted to pseudobulk, logCPM, TMM-normalized expression, or z-scores. The full available gene set is retained rather than only the Youth Score signature genes.

This .h5ad file is intended for Geneformer or other cell-level models that require raw counts. It is different from the pseudobulk logCPM input used by the Youth Score parser.

Important: The data are from Mus musculus. Before running Geneformer, confirm whether the selected checkpoint supports mouse genes directly or requires mouse-to-human ortholog mapping.

The training cohort has only two young mice, both female. Therefore, age and sex are partially confounded. The Youth Score feature-selection pipeline adjusts for sex, but neither the raw dataset nor the final score should be described as sex-independent.
