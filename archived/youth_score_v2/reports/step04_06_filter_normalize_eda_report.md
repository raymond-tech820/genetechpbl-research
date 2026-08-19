# FACS v2 Step 04-06: Gene Filtering, TMM Normalization, and EDA

## Scope

This full-data normalization is for QC and EDA only. Nested validation and model training must repeat filtering and normalization inside each training fold.

No files under `data_facs` were modified.

## Gene Filter

- Raw genes: 22966
- Retained genes: 12859
- Removed genes: 10107
- Rule: CPM > 1 in at least 2 mice

## TMM

- Samples: 14
- Median raw library size: 94384073
- Median TMM effective library size: 87128862.87

## PCA

- PCA used top variable genes: 2000
- PCA did not use all retained genes.

## PC Associations

  pc spearman_age_months spearman_sex_male_numeric spearman_n_cells
 PC1           0.3469174                0.54912518       -0.4268429
 PC2          -0.6516421                0.07844645       -0.5720576
 PC3          -0.3422293                0.19611614       -0.2574259
 spearman_raw_library_size spearman_effective_library_size_tmm
               -0.09010989                         -0.39780220
               -0.78461538                         -0.63956044
               -0.21758242                         -0.05494505
 spearman_detected_genes
             -0.98241758
              0.08571429
              0.24395604

## Outputs

- `outputs/facs_v2/processed/facs_v2_limb_msc_full_data_filtered_counts.rds`
- `outputs/facs_v2/processed/facs_v2_limb_msc_full_data_tmm_logcpm.rds`
- `outputs/facs_v2/processed/facs_v2_limb_msc_full_data_dge_tmm.rds`
- `outputs/facs_v2/processed/facs_v2_limb_msc_full_data_voom.rds`
- `outputs/facs_v2/qc/step04_full_data_gene_filter_qc.csv`
- `outputs/facs_v2/qc/step05_tmm_normalization_factors.csv`
- `outputs/facs_v2/eda/step06_pca_scores_high_variable_genes.csv`
- `outputs/facs_v2/eda/step06_pc_associations.csv`
- `outputs/facs_v2/eda/step06_pca_by_age_group.png`
- `outputs/facs_v2/eda/step06_pca_by_sex.png`
- `outputs/facs_v2/eda/step06_pca_by_raw_library_size.png`
