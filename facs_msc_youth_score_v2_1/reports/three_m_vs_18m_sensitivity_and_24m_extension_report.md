# FACS v2.1 Step 15: 3m-vs-18m Structural Sensitivity and 24m Extension

## Scope

This step removes the male-only 24m group from training to test whether pooled-old design drives instability. The frozen 3m/18m model is then applied to 24m male mice without refitting or recalibration.

No files under `data_facs` were modified.

## 3m-vs-18m Nested LOMO Summary

 n_mice       auc    all_rho library_rho detected_rho cell_count_rho
     10 0.8333333 -0.5685352  -0.3333333    0.3818182      0.1030303
 young_minus_old_median zero_signature_folds full_rank_folds weak_support_folds
              0.1190158                    0              10                  6
 median_jaccard
      0.3333333

## Effect Correlation vs Pooled 18/24m Model

 common_effect_spearman female_effect_spearman male_effect_spearman
              0.8620589              0.9999565            0.7883142
 compared_genes
          12411

## Frozen 3m/18m -> 24m Male Extension

 age age_months    sex       score
  3m          3 female  1.03105607
 18m         18 female  0.01791279
  3m          3   male  0.96626106
 18m         18   male -0.09907982
 24m         24   male  0.42650410

## Male Age Trajectory Check

 male_age_spearman median_3m_male median_18m_male median_24m_male
        -0.6227992      0.9662611     -0.09907982       0.4265041
 expected_3_gt_18_gt_24
                  FALSE

## Outputs

- `three_m_vs_18m_nested_lomo_scores.csv`
- `three_m_vs_18m_fold_signatures.csv`
- `three_m_vs_18m_gene_effects_by_fold.csv`
- `three_m_vs_18m_summary.csv`
- `three18_frozen_model_24m_extension_scores.csv`
- `three18_frozen_model_24m_extension_male_summary.csv`
- `three_m_vs_18m_effect_correlation_vs_pooled.csv`
- `three_m_vs_18m_signature_overlap_vs_pooled.csv`
