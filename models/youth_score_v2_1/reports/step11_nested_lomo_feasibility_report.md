# FACS v2.1 Step 11: Fully Nested LOMO Feasibility Audit

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
       -0.7500916             0.2182179     1.0000000    0.5960937  0.2822987
       -0.6328898             0.6546537     0.9791667    0.5412196  0.2781036
       -0.7125870             0.4364358     1.0000000    0.5487876  0.2659956
 score_cor_library score_cor_detected score_cor_cell_count median_gene_coverage
        -0.6527473        -0.03736264           -0.5544558                    1
        -0.5208791        -0.11648352           -0.4620465                    1
        -0.6439560         0.17362637           -0.4070410                    1
 median_weighted_gene_coverage median_ribosomal_fraction
                             1                         0
                             1                         0
                             1                         0
 median_cell_cycle_fraction median_gene_library_burden
                          0                  0.6062024
                          0                  0.5719833
                          0                  0.5312285
 median_gene_detected_burden median_gene_cell_count_burden
                   0.2022009                     0.4169876
                   0.2296322                     0.3692910
                   0.2432032                     0.3405693

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

- `outputs_final/facs_v2_1/validation/step11_nested_lomo_scores.csv`
- `outputs_final/facs_v2_1/validation/step11_nested_lomo_fold_diagnostics.csv`
- `outputs_final/facs_v2_1/validation/step11_nested_lomo_fold_signatures.csv`
- `outputs_final/facs_v2_1/validation/step11_nested_lomo_gene_effects_by_fold.csv`
- `outputs_final/facs_v2_1/validation/step11_nested_lomo_fold_pc_technical_audit.csv`
- `outputs_final/facs_v2_1/validation/step11_nested_lomo_signature_jaccard.csv`
- `outputs_final/facs_v2_1/validation/step11_nested_lomo_gene_selection_frequency.csv`
- `outputs_final/facs_v2_1/validation/step11_lowest_depth_removal_sensitivity.csv`
