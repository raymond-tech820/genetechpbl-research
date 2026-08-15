# FACS Youth Score v2.1 Step 17: Formal Nested Permutation Null for Primary Medium Model

Attempted permutations requested: 999
Valid null permutations completed: 999
Failed permutations: 0

## Locked Definition

- Model: Factorial Medium Original only.
- Permutation: complete 3m/18m/24m age labels are reassigned within sex strata, then `age_months` and `age_group` are regenerated from the permuted full age label.
- Female label composition is preserved as 2 x 3m and 2 x 18m.
- Male label composition is preserved as 4 x 3m, 2 x 18m, and 4 x 24m.
- Every permutation reruns the full nested LOMO training-side workflow, including gene filtering, TMM, voom, factorial DE, ranking, signature selection, weights, calibration, and held-out scoring.
- Primary statistic: `abs_all_age_rho = abs(Spearman(score, age_months))`.
- Supporting statistic: `auc_distance_from_half = abs(AUC - 0.5)`.
- Empirical p-value: `(1 + number of null statistics >= observed statistic) / (valid_null_n + 1)`.

This is a practical training-pipeline permutation null. It tests whether the predefined nested factorial pipeline produces stronger age-label association than sex-stratified randomized age labels. It does not test technical independence, external validity, or whether the score is a continuous biological aging clock.

## Summary

| Metric | Role | Observed | Null median | Null 95% interval | Empirical p | Valid null n |
|---|---|---:|---:|---:|---:|---:|
| abs_all_age_rho | primary | 0.633 | 0.352 | [0.019, 0.811] | 0.153 | 999 |
| auc_distance_from_half | supporting | 0.479 | 0.208 | [0.000, 0.479] | 0.039 | 999 |
| young_minus_old_median | recorded_not_primary | 0.263 | -0.081 | [-0.277, 0.153] | 0.002 | 999 |

## Diagnostics

- Total metric rows: 1000
- Duplicate assignment hashes among valid null permutations: 0
- Observed zero-signature folds: 0
- Observed rank-deficient folds: 0

## Outputs

- `nested_permutation_medium_scores.csv`
- `nested_permutation_medium_fold_diagnostics.csv`
- `nested_permutation_medium_metrics.csv`
- `nested_permutation_medium_summary.csv`
