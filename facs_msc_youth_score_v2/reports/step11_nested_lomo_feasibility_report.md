# FACS v2 Step 11: Fully Nested LOMO Feasibility Audit

## Scope

This is the first fully nested LOMO run for feasibility auditing. It is not the final model-selection or permutation stage.

Every outer fold repeats training-side gene filtering, TMM, voom, age-only DE, factorial DE, contrasts, sex-linked exclusion, reliability audit, ranking, signature selection, weights, frozen deployable normalization, and held-out scoring.

No files under `data_facs` were modified.

## Main Feasibility Summary

 signature_version n_folds zero_signature_folds median_signature_size
             Small      14                    0                    40
            Medium      14                    0                   100
             Large      14                    0                   200
 all_spearman_age old_only_spearman_age young_old_auc young_median old_median
       -0.7688439             0.1091089             1    0.6376526  0.3013938
       -0.7125870             0.4364358             1    0.5348597  0.2700303
       -0.6938347             0.5455447             1    0.5247365  0.3052094
 score_cor_library score_cor_detected score_cor_cell_count median_gene_coverage
        -0.8241758        -0.12967033           -0.7678773                    1
        -0.8461538        -0.04175824           -0.7106715                    1
        -0.8417582        -0.03296703           -0.6754680                    1
 median_weighted_gene_coverage median_ribosomal_fraction
                             1                     0.000
                             1                     0.000
                             1                     0.005
 median_cell_cycle_fraction median_gene_library_burden
                          0                  0.6383908
                          0                  0.5912563
                          0                  0.5593503
 median_gene_detected_burden median_gene_cell_count_burden
                   0.2146131                     0.4487836
                   0.2458054                     0.3964656
                   0.2637817                     0.3772293

## Fold Diagnostics

- Folds: 14
- Full-rank factorial folds: 14 / 14
- Weak-support folds: 4 / 14
- Zero-signature folds across all versions: 0
- Lowest-depth mouse for sensitivity: 3_9_M

## Technical Axis Audit

Per-fold PC1/PC2 technical correlations are recorded in `step11_nested_lomo_fold_pc_technical_audit.csv`.
Per-score held-out technical associations are summarized above and retained per mouse in `step11_nested_lomo_scores.csv`.

## Weak-Support Fold Handling

Weak-support folds are not removed. They are flagged for separate interpretation because one age-by-sex training cell has only one mouse.

## Stability Gate Candidates

This run records gene-level fold effects and selection frequency. Final primary candidates should require high common-direction consistency, male direction stability, no systematic female reversal, age-only/common agreement, and adequate fold-selection frequency.

## Outputs

- `outputs/facs_v2/validation/step11_nested_lomo_scores.csv`
- `outputs/facs_v2/validation/step11_nested_lomo_fold_diagnostics.csv`
- `outputs/facs_v2/validation/step11_nested_lomo_fold_signatures.csv`
- `outputs/facs_v2/validation/step11_nested_lomo_gene_effects_by_fold.csv`
- `outputs/facs_v2/validation/step11_nested_lomo_fold_pc_technical_audit.csv`
- `outputs/facs_v2/validation/step11_nested_lomo_signature_jaccard.csv`
- `outputs/facs_v2/validation/step11_nested_lomo_gene_selection_frequency.csv`
- `outputs/facs_v2/validation/step11_lowest_depth_removal_sensitivity.csv`
