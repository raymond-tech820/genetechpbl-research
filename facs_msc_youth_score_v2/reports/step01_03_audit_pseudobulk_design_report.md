# FACS v2 Step 01-03: Input Audit, Pseudobulk, and Design Audit

## Inputs

- Metadata: `/Users/strangenoah/Desktop/AI+X/Genetech/Youth_score/data_facs/limb_muscle_msc/facs_limb_muscle_msc_young_old_metadata.csv`
- Matrix: `/Users/strangenoah/Desktop/AI+X/Genetech/Youth_score/data_facs/limb_muscle_msc/expression_bpcells_young_old`

No files under `data_facs` were modified by this script. All new outputs were written under `outputs/facs_v2`.

## Input Read Check

- Matrix dimensions: 22966 genes x 935 cells
- Metadata rows: 935
- Matrix columns match metadata `index`: TRUE
- Metadata `n_counts` exact matches: 933 / 935
- Metadata `n_genes` exact matches: 935 / 935

## Cohort Structure

   mouse age age_group    sex n_cells pseudobulk_library_size
  3_38_F  3m     Young female      77                79934534
  3_39_F  3m     Young female      53                66348808
  3_10_M  3m     Young   male      28                28380229
  3_11_M  3m     Young   male      14                17415067
   3_8_M  3m     Young   male      75                82136470
   3_9_M  3m     Young   male      17                 9599638
 18_46_F 18m       Old female     157               261459669
 18_47_F 18m       Old female     121               296683550
 18_45_M 18m       Old   male      31               108019124
 18_53_M 18m       Old   male      45                87537694
 24_58_M 24m       Old   male      79               132161042
 24_59_M 24m       Old   male      68               161460489
 24_60_M 24m       Old   male      91               116295606
 24_61_M 24m       Old   male      79               101230452
 pseudobulk_detected_genes
                     12876
                     12415
                     12456
                     10153
                     14241
                     11393
                     14857
                     14205
                     11115
                     12432
                     11177
                     10615
                     11216
                     11634

## Age Summary

 age_group age age_months n_cells n_mice
     Young  3m          3     264      6
       Old 18m         18     354      4
       Old 24m         24     317      4

## Full-Data Design Diagnostics

                                design n_samples n_columns rank full_rank
                   age_only_~age_group        14         2    2      TRUE
               additive_~sex+age_group        14         3    3      TRUE
 factorial_cell_means_~0+sex_age_group        14         4    4      TRUE
 residual_df condition_number
          12         2.820934
          11         4.042026
          10         1.732051
                                                                                         columns
                                                                        (Intercept);age_groupOld
                                                                (Intercept);sexmale;age_groupOld
 age_sex_groupfemale_Young;age_sex_groupmale_Young;age_sex_groupfemale_Old;age_sex_groupmale_Old

## Fold-Level Design Notes

- Outer LOMO folds audited: 14
- Full-rank factorial folds: 14 / 14
- Folds with any age-by-sex cell supported by <2 training mice: 4 / 14

Folds holding out one of the two female Young or female Old mice remain algebraically estimable but are weakly supported for sex-specific effects.

## Outputs

- `outputs/facs_v2/processed/facs_v2_limb_msc_mouse_metadata.csv`
- `outputs/facs_v2/processed/facs_v2_limb_msc_pseudobulk_counts.csv`
- `outputs/facs_v2/processed/facs_v2_limb_msc_pseudobulk_counts.rds`
- `outputs/facs_v2/qc/step01_cell_qc_matrix_metadata_check.csv`
- `outputs/facs_v2/qc/step01_age_summary.csv`
- `outputs/facs_v2/qc/step01_age_by_sex_mouse_counts.csv`
- `outputs/facs_v2/qc/step01_age_by_sex_cell_counts.csv`
- `outputs/facs_v2/qc/step03_full_design_diagnostics.csv`
- `outputs/facs_v2/qc/step03_fold_design_diagnostics.csv`
