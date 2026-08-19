# FACS v2.1 Step 01-03: Raw Data Audit, Pseudobulk, and Design Audit

## Scope

The FACS v2.1 model starts from the diaphragm-excluded `data_facs/raw_data` bundle. The raw files are read only; all outputs are written under `outputs_final/facs_v2_1`.

## Raw Data Cleanliness

- Matrix dimensions: 22966 genes x 815 cells
- Metadata rows: 815
- Mice: 14
- Young 3m cells: 264
- Old 18m/24m cells: 551
- Tissue == Limb_Muscle: 815 / 815
- MSC ontology: 815 / 815
- Diaphragm subtissue matches: 0
- Standardized fore/hind limb subtissue cells: 815 / 815

## Original Subtissue Values


       forelimb and hindlimb          ForelimbandHindlimb 
                          31                          581 
Muscle forelimb and hindlimb 
                         203 

## Donor Pseudobulk Summary

   mouse age age_group    sex         subtissue n_cells pseudobulk_library_size
  3_38_F  3m     Young female Forelimb_Hindlimb      77                79934534
  3_39_F  3m     Young female Forelimb_Hindlimb      53                66348808
  3_10_M  3m     Young   male Forelimb_Hindlimb      28                28380229
  3_11_M  3m     Young   male Forelimb_Hindlimb      14                17415067
   3_8_M  3m     Young   male Forelimb_Hindlimb      75                82136470
   3_9_M  3m     Young   male Forelimb_Hindlimb      17                 9599638
 18_46_F 18m       Old female Forelimb_Hindlimb      82                97642498
 18_47_F 18m       Old female Forelimb_Hindlimb      76                89804246
 18_45_M 18m       Old   male Forelimb_Hindlimb      31               108019124
 18_53_M 18m       Old   male Forelimb_Hindlimb      45                87537694
 24_58_M 24m       Old   male Forelimb_Hindlimb      79               132161042
 24_59_M 24m       Old   male Forelimb_Hindlimb      68               161460489
 24_60_M 24m       Old   male Forelimb_Hindlimb      91               116295606
 24_61_M 24m       Old   male Forelimb_Hindlimb      79               101230452
 pseudobulk_detected_genes
                     12876
                     12415
                     12456
                     10153
                     14241
                     11393
                     13280
                     13216
                     11115
                     12432
                     11177
                     10615
                     11216
                     11634

## Design Audit

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

## Fold Support

Weak-support folds are retained and flagged rather than silently changing the model.

 heldout_mouse heldout_age heldout_sex train_n_female_Young train_n_male_Young
        3_38_F          3m      female                    1                  4
        3_39_F          3m      female                    1                  4
        3_10_M          3m        male                    2                  3
        3_11_M          3m        male                    2                  3
         3_8_M          3m        male                    2                  3
         3_9_M          3m        male                    2                  3
       18_46_F         18m      female                    2                  4
       18_47_F         18m      female                    2                  4
       18_45_M         18m        male                    2                  4
       18_53_M         18m        male                    2                  4
       24_58_M         24m        male                    2                  4
       24_59_M         24m        male                    2                  4
       24_60_M         24m        male                    2                  4
       24_61_M         24m        male                    2                  4
 train_n_female_Old train_n_male_Old design_full_rank
                  2                6             TRUE
                  2                6             TRUE
                  2                6             TRUE
                  2                6             TRUE
                  2                6             TRUE
                  2                6             TRUE
                  1                6             TRUE
                  1                6             TRUE
                  2                5             TRUE
                  2                5             TRUE
                  2                5             TRUE
                  2                5             TRUE
                  2                5             TRUE
                  2                5             TRUE
 weakly_supported_factorial_cell
                            TRUE
                            TRUE
                           FALSE
                           FALSE
                           FALSE
                           FALSE
                            TRUE
                            TRUE
                           FALSE
                           FALSE
                           FALSE
                           FALSE
                           FALSE
                           FALSE

## Outputs

- `outputs_final/facs_v2_1/processed/facs_v2_1_limb_msc_mouse_metadata.csv`
- `outputs_final/facs_v2_1/processed/facs_v2_1_limb_msc_pseudobulk_counts.rds`
- `outputs_final/facs_v2_1/qc/step01_03_raw_data_validation.csv`
- `outputs_final/facs_v2_1/qc/step03_fold_design_diagnostics.csv`
