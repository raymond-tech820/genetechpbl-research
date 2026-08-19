# FACS v2.1 Step 04-06: Gene Filtering, TMM Normalization, and EDA

## Scope

This full-data normalization is for QC and EDA only. Nested validation and model training must repeat filtering and normalization inside each training fold.

No files under `data_facs` were modified.

## Gene Filter

- Raw genes: 22966
- Retained genes: 12824
- Removed genes: 10142
- Rule: CPM > 1 in at least 2 mice

## TMM

- Samples: 14
- Median raw library size: 88670970
- Median TMM effective library size: 89248784.96

## PCA

- PCA used top variable genes: 2000
- PCA did not use all retained genes.

## PC Associations

  pc spearman_age_months spearman_sex_male_numeric spearman_n_cells
 PC1           0.4172384                0.47067872       -0.1936195
 PC2          -0.6516421                0.07844645       -0.5984602
 PC3          -0.4313027                0.15689291       -0.3190321
 spearman_raw_library_size spearman_effective_library_size_tmm
                 0.3406593                         -0.22637363
                -0.8329670                         -0.66153846
                -0.3846154                         -0.09450549
 spearman_detected_genes
             -0.95164835
              0.09450549
              0.20000000

## Outputs

- `outputs_final/facs_v2_1/processed/facs_v2_1_limb_msc_full_data_filtered_counts.rds`
- `outputs_final/facs_v2_1/processed/facs_v2_1_limb_msc_full_data_tmm_logcpm.rds`
- `outputs_final/facs_v2_1/processed/facs_v2_1_limb_msc_full_data_dge_tmm.rds`
- `outputs_final/facs_v2_1/processed/facs_v2_1_limb_msc_full_data_voom.rds`
- `outputs_final/facs_v2_1/qc/step04_full_data_gene_filter_qc.csv`
- `outputs_final/facs_v2_1/qc/step05_tmm_normalization_factors.csv`
- `outputs_final/facs_v2_1/eda/step06_pca_scores_high_variable_genes.csv`
- `outputs_final/facs_v2_1/eda/step06_pc_associations.csv`
- `outputs_final/facs_v2_1/eda/step06_pca_by_age_group.png`
- `outputs_final/facs_v2_1/eda/step06_pca_by_sex.png`
- `outputs_final/facs_v2_1/eda/step06_pca_by_raw_library_size.png`
