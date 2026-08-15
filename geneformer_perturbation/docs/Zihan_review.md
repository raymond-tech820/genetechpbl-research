# Zihan: `limb_msc_youth_score_v1` (TMS Droplet Limb Muscle MSC, pseudobulk gene signature) — Review

## 1. What was built

- **Data.** 12 mice pseudobulked from TMS Droplet Limb Muscle MSCs (9,649 cells). Young = 3 mo (2 mice, **both female**); Old = 18/21/24 mo (10 mice, 4 F + 6 M). 12,735 genes after CPM>1-in-≥2 filter, TMM normalized (edgeR).
- **DE.** edgeR quasi-likelihood, design `~ sex + age_group` on all 12 mice.
- **Ranking.** `q_g = |adjusted_logFC| · |age_rho| · r_g` with `r_g = π_LOMO · π_depth · π_sex` (π_depth and π_sex are 0/1 gates, π_LOMO is the LOMO sign-agreement rate).
- **Signatures.** Small (20+20), Medium (50+50, selected as primary), Large (100+100, comparator). Score = weighted mean z(young-module) − weighted mean z(old-module), calibrated by training medians so young median = 1, old median = 0.
- **Validation.** Nested LOMO — every step (gene filter, TMM, DE, LOMO stability, low-depth check, sex filter, ranking, gene selection, mean/SD, calibration) is refit inside each fold.
- **Robustness.** 20 mouse-level age-label permutations; 300 expression-matched random gene sets; 1000 weight shuffles; low-depth (`18-F-50`, `18-F-51`) drop test; no-sex-filter stress; single-cell aggregate vs pseudobulk (Spearman 0.993).

## 2. What is genuinely strong

1. **Leakage control is textbook-clean.** The nested-LOMO code (`validate_limb_muscle_msc_step14_nested_lomo.R`) refits *everything* per fold — including TMM using a training-only reference sample — with no full-data leakage into the held-out score. Very few undergraduate/PBL projects get this right.
2. **Documentation of limitations is exemplary.** The report explicitly labels young-heldout folds as `held_out_young_stress_test`, `bootstrap_stability = not_assessed`, `permutation_scope = practical_training_pipeline_null`. These flags are also encoded in `calibration.json` — the model file itself carries its own caveats. This is best-practice.
3. **Multi-layer control battery.** Age-label permutation + expression-matched gene set + weight shuffle + low-depth drop + no-sex-filter stress covers most of the "is this signal real?" attack surface. The finding that weight shuffling barely changes performance (i.e. the signal lives in *gene identity + module direction*, not the fine-grained weights) is the correct read.
4. **Coverage QC in the parser.** `score_limb_msc_youth.R` gates on `gene_coverage` and `weighted_coverage` against the calibration JSON — a scorer that fails safely on partial vocabularies.
5. **Parser numerical equivalence to Step 12** at ~1e-15 — the shipped R function reproduces the training scores bit-for-bit.

## 3. What is fragile — in decreasing order of impact

### (a) The young-female confound is structural, not just documentational

Two young mice, both female, versus ten old mice of mixed sex means `~ sex + age_group` can estimate a sex main effect *only from old males vs old females* and an age effect *only from young-female vs (sex-adjusted old)*. Any **sex × age interaction** — biologically very likely for MSC populations — is aliased into the age coefficient and cannot be diagnosed. The `π_sex` filter (drop genes where |old-male − old-female| > 1 and > 2× |age logFC|) helps only for genes that *already look sex-sensitive in the old group*; a gene that is age-driven only in females passes silently. The report acknowledges "not sex-independent"; the report should go further and say the age effect is estimated **within young females** and generalisation to young males is unmodelled.

### (b) The headline ρ = −0.94 is driven by old-fold ordering, not young–old separation

From `step14_nested_lomo_summary.csv` (Medium):

| slice | rows | Spearman(score, age) |
|---|---:|---:|
| all mice | 12 | **−0.939** |
| old-only | 10 | **−0.934** |

The old-only correlation is essentially identical to the all-mice one — the model is largely ranking 18 vs 21 vs 24 mo old mice among themselves, then adding two well-separated young points on top. In `step15_lomo_score_by_age.png` the two red (young stress-test) points sit at 0.35–0.50, but the top 18 mo old-male triangles reach 0.30–0.44 (circle = female, triangle = male in the plot legend) — young–old separation is real but the overlap is not trivial. Report should present **old-only ρ as the primary evidence** and treat the whole-cohort ρ as an artifact of the design.

### (c) The permutation p-value is at its floor

Only 20 permutations, so the minimum attainable empirical p is (0+1)/(20+1) ≈ 0.048. The observed |ρ| = 0.94 sits far above the null in `step16_age_label_permutation_null.png`, but the reported p-value is architectural (it *can't be smaller*), not evidentiary. **Rerun with ≥ 999 permutations** (cheap — a few CPU-hours at most) so the p-value is actually informative. Better still, do a **nested outer-LOMO permutation** (permute mouse-level labels, then rerun the full outer LOMO) rather than the "practical training-pipeline null" — that is what the report itself flags as future work.

### (d) `bootstrap_stability = not_assessed` should be closed

Per-mouse cell counts range from 271 to 1,328 (median ~840). A **cell-level bootstrap within each mouse → repseudobulk → rescore** loop would produce a mouse-level score CI at essentially zero cost and would let the deliverable put an uncertainty bar on every mouse in Fig. `step15_lomo_score_by_age.png`. This is the single most valuable easy addition.

### (e) TMM held-out normalization has a small leakage seam

`heldout_tmm_logcpm()` builds a combined DGE with training + one held-out sample and refits TMM. Even with a fixed training-median `refColumn`, the effective library sizes of the *training* samples shift slightly when the held-out sample is added. Cleaner: compute training-only TMM factors, then normalise the held-out mouse using edgeR's `getCounts` + fixed offset (or `edgeR::processAmplicons`-style external offset). Effect on ρ is likely tiny, but it's a defensible cleanup.

### (f) Weight formula is ad hoc

`w_g = min(|adjusted_logFC|, 3) · r_g`. Because π_depth and π_sex are 0/1, `r_g` collapses to π_LOMO for surviving genes, so `w_g ≈ min(|logFC|, 3) · π_LOMO`. The clip at 3 has no obvious rationale beyond "big effects shouldn't dominate"; **document the choice or replace with a shrinkage estimator** (empirical-Bayes moderated logFC from `limma::eBayes` is a principled equivalent).

### (g) Calibration is training-median rescaling, not probability calibration

By construction training-set young median = 1 and old median = 0. Held-out young scores land around 0.4 (see step15 figure), so the [0,1] scale is not a probability and not preserved on unseen data. **Rename `youth_score_clipped_0_1` to `youth_score_scaled_relative_to_training_medians`** or similar; alternatively fit a Platt/isotonic calibrator on the LOMO out-of-fold scores against Young/Old labels so downstream users get an actual probability.

### (h) Signature stability is moderate

Median cross-fold Jaccard is 0.54 for Medium (0.58 for Large). About half of the 100 genes turn over between LOMO folds — not disastrous, but the "final Medium signature" is one of many nearly-equivalent choices. Consider reporting **stability-selected genes** (present in ≥ k of 12 folds) as a secondary, more defensible list.

### (i) Missing baselines

No comparison to (i) first PC of TMM logCPM, (ii) a published mouse aging clock (e.g. Meer 2018, Petkovich 2017), (iii) a naive elastic-net regression on all filtered genes. Without at least a 1-PC baseline the ρ = −0.94 has no reference for how impressive it is on a 12-mouse cohort. **Add PC1 and elastic-net baselines to Section 3.** Kaile's SCAT elastic-net Youth Score already exists — that's your baseline.

## 4. Smaller notes

- Report writes `age_de$logFC` as "log2FC" in the volcano axis but as `logFC` elsewhere. edgeR returns log2FC by default — clarify once, uniformly.
- `step16_age_label_permutation_null.png` uses a legend "True Large/Medium" with vertical lines that overlap near |ρ| ≈ 0.94; the two truth lines are indistinguishable. Add horizontal labels or split into two facets.
- `pmax(abs(...), 1e-6)` guards in the sex-sensitivity ratio are fine but should be documented alongside the threshold `2` — right now the ratio for near-zero logFC genes goes to ∞ and is trimmed silently.
- `check.names = FALSE` is used consistently — good.
- `set.seed(20260717)` only in Step 16; consider setting seeds everywhere randomness enters (permutation, random gene sets, weight shuffles).

## 5. Recommended next steps

1. **Rerun age-label permutation with ≥ 999 iterations**, nested-outer-LOMO if compute allows.
2. **Add cell-level bootstrap** per mouse to close `bootstrap_stability`.
3. **Report old-only ρ as the primary metric** and demote whole-cohort ρ to a supplementary panel.
4. **Add PC1 and elastic-net baselines** (borrow Kaile's).
5. **Recruit young males if possible.** No amount of statistics fixes the young-female confound; only new samples do. If TMS has any young-male Droplet data left, use it — even one extra mouse would be worth more than any further robustness control.
6. **Cross-check against Kaile's frozen Limb-FACS model** on the same 12 Droplet mice — if two independently trained signatures agree at the mouse level, that is stronger evidence than either alone.
