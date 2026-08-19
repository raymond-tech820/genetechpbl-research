# FACS v2 Step 07-10: Full-Data DE, Sex Consistency, and Candidate Ranking

## Scope

This is a full-data exploratory analysis used to audit signal structure and generate provisional candidate rankings. It is not held-out validation evidence.

Nested LOMO model training must repeat gene filtering, TMM, voom, DE, sex-consistency filtering, and ranking inside each training fold.

No files under `data_facs` were modified.

## Models

- Age-only baseline DE: `~ age_group`.
- Primary exploratory DE: `~ 0 + sex:age_group` cell-means factorial design.
- Female age effect: female Old - female Young.
- Male age effect: male Old - male Young.
- Common age effect: average of female and male age effects.
- Interaction: male age effect - female age effect.

All logFC fields are log2 fold-changes with old - young orientation.

## Candidate Reliability Rules

- Female and male age effects must have the same nonzero direction.
- Age-only and common effects must have the same nonzero direction.
- Hard-coded obvious sex-linked genes are excluded.
- Interaction magnitude uses epsilon = 1e-06.
- Ranking score: `abs(common_log2FC) * abs(common moderated t) * interaction_penalty`.

## Summary

                             metric value
              filtered_genes_tested 12859
                age_only_fdr_lt_0_1   287
                  common_fdr_lt_0_1    77
                  female_fdr_lt_0_1    11
                    male_fdr_lt_0_1   623
             interaction_fdr_lt_0_1     1
     sex_direction_concordant_genes  7457
   age_only_common_concordant_genes 11995
 primary_reliability_pass_full_data  7443
       young_high_ranked_candidates  4034
         old_high_ranked_candidates  3409
           sex_linked_genes_present     7

## Interpretation Notes

- FDR-significant genes are not the sole definition of the final signature.
- Candidate ranking uses effect size, moderated statistic, sex-direction concordance, age-only/common agreement, and interaction penalty.
- Stability, low-depth sensitivity, and nested held-out performance are not assessed in this step.

## Outputs

- `outputs/facs_v2/de/step07_10_full_data_factorial_de_and_ranking.csv`
- `outputs/facs_v2/de/step07_10_full_data_de_summary.csv`
- `outputs/facs_v2/signature/step10_full_data_young_high_ranked_candidates.csv`
- `outputs/facs_v2/signature/step10_full_data_old_high_ranked_candidates.csv`
- `outputs/facs_v2/signature/step10_full_data_top100_each_direction_candidates.csv`
- `outputs/facs_v2/de/step07_common_vs_interaction_effects.png`
- `outputs/facs_v2/de/step07_female_vs_male_age_effects.png`
