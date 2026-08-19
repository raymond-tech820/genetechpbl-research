# Youth Score v3.1 M1-M4 External Comparison: Frozen Pipeline

Status: frozen before M1/M4 external scoring; amended after input-unit provenance audit

Freeze date: 2026-08-05

This document freezes the comparison of the FACS v3.1 M1 and M4 models. It does not authorize model redevelopment, threshold tuning, or outcome-driven feature changes.

## 1. Scientific Questions

The comparison asks three distinct questions:

1. Training-cohort statistical evidence: does each fully nested pipeline produce stronger age-label association than its sex-stratified randomized-label null?
2. Biological external application: do the frozen scores produce interpretable baseline and treatment associations in GSE176206?
3. Cross-assay transportability: do the frozen FACS-derived scores preserve direction and donor ordering in TMS Droplet Limb Muscle MSC?

These questions must not be collapsed into a single claim of model validity.

## 2. Frozen Models and Roles

### M1

Frozen ID: `M1_raw_all_gene_denominator`

Role: primary deployable reference.

- `factorial_stability_selected` framework;
- inner selection frequency threshold $\pi_g\geq0.75$;
- 58 genes: 30 young-high and 28 old-high;
- raw all-gene $\log_2(\mathrm{CPM}+1)$ deployment transform;
- frozen FACS gene means, standard deviations, weights, modules, and calibration centers.

### M4

Frozen ID: `M4_raw_all_gene_denominator_pi_0_90`

Role: secondary, post-ablation exploratory candidate.

- identical to M1 except $\pi_g\geq0.90$;
- 29 genes: 16 young-high and 13 old-high;
- M4 is a strict subset of M1 in the full-data export.

M4 was selected after examining the controlled ablations. A favorable external result may strengthen its candidate role but cannot erase this selection history.

### M2 and M3

M2 and M3 are excluded from the primary frozen external comparison. M2 remains an archived negative result for the simple DE technical-adjustment hypothesis; M3 remains an archived threshold sensitivity. Neither may be added after viewing M1/M4 external outcomes.

A future M2 deployment analysis requires a separately frozen, explicitly named parser and a separate protocol amendment made before inspecting its results.

## 3. Completed Internal Evidence

The following steps are complete and are not rerun or redefined:

- full-data M1/M4 artifact export;
- parser numerical equivalence and coverage-gate tests;
- fully nested observed LOMO reconstruction;
- 20-permutation hardened smoke test;
- 999 unique sex-stratified complete-age permutations.

Formal primary results:

| Model | Observed $|\rho_{all-age}|$ | Null mean | $Z_{null}$ | Empirical $p$ |
|---|---:|---:|---:|---:|
| M1 | 0.8063 | 0.3585 | 1.953 | 0.023 |
| M4 | 0.8439 | 0.3547 | 2.177 | 0.010 |

These results support non-random age-label association under the frozen practical null. They do not establish technical independence, external validity, a continuous aging clock, or M4 superiority.

## 4. Shared Frozen Inference Rules

Both external datasets must use the same model artifacts and parser implementation.

The input path is:

$$
\text{single-cell raw counts}
\rightarrow
\text{biological-unit pseudobulk}
\rightarrow
\text{frozen raw-library transform}
\rightarrow
\text{frozen gene standardization}
\rightarrow
\text{frozen module score and calibration}.
$$

For gene $g$ and biological unit $j$:

$$
x_{gj}=\log_2\left(10^6\frac{C_{gj}}{L_j}+1\right),
\qquad
z_{gj}=\frac{x_{gj}-\mu_g^{FACS}}{s_g^{FACS}}.
$$

The raw score is:

$$
S_j=
\frac{\sum_{g\in G_Y}w_gz_{gj}}{\sum_{g\in G_Y}w_g}
-
\frac{\sum_{g\in G_O}w_gz_{gj}}{\sum_{g\in G_O}w_g}.
$$

The calibrated score is:

$$
Y_j=\frac{S_j-M_O^{FACS}}{M_Y^{FACS}-M_O^{FACS}}.
$$

External data must not alter genes, weights, modules, $\mu_g$, $s_g$, calibration centers, orientation, thresholds, or missing-gene rules.

## 5. Coverage Gate

Coverage is audited before outcome metrics. Report separately for each model:

- detected and missing signature genes;
- unweighted and weighted total coverage;
- young-high unweighted and weighted coverage;
- old-high unweighted and weighted coverage;
- effective usable signature size.

The frozen minimum is 0.80 for total weighted coverage and for each module's weighted coverage. A model below any threshold is labeled `feature_limited`; its biological or transportability performance is not interpreted.

M4's smaller signature makes this check particularly important. No missing gene may be replaced.

## 6. Phase 1: GSE176206 Biological External Application

GSE176206 is run first because it is closest to the regenerative-perturbation use case. It is an independent-dataset biological application and directional validation, not proof of a universal Youth Score or chronological aging clock.

### 6.1 Input and units

Use `data_geo/processed/GSE176206_msc_sokm.h5ad` as the authoritative expression source and use only `layers['counts']`. The frozen source contract is 20,661 cells by 28,694 genes with the required `age`, `treatment`, `animal`, `sample`, and `state` metadata. The files under `data_geo/processed/GSE176206_parsed` are program-friendly derivatives and may be used only after counts, genes, cells, and grouping equivalence are verified against that H5AD.

Both models use the same cells, gene mapping, duplicate-gene handling, metadata exclusions, and pseudobulk matrix.

Primary units are known biological animals nested within age and exact treatment arm:

$$
C_{g,a,t,i}=\sum_{c\in(a,t,i)}C_{gc}.
$$

The exact treatment arms are `Tg+/Dox+`, `Tg+/Dox-`, and `Tg-/Dox+`. `Tg+/Dox+` is the SOKM arm. The other two arms are distinct controls and must not be merged by matching their numeric `animal` labels.

The current data support 18 known `age x exact-treatment x animal` units: two ages, three treatment arms, and three known animals per arm. Cells from multiple `sample` or library rows may be summed only when age, exact treatment, and nested animal label are all identical.

The original study states that animal labels are specific to each `age:treatment` combination because those combinations were run in separate library-preparation reactions. Therefore equal numeric labels across treatment arms do not identify the same animal. Control-versus-SOKM comparisons are unpaired. Rows with `animal = unknown` are excluded from primary summaries. This unit definition is grounded in the GSE176206 study design and the Methods of Roux et al., Cell Systems 2022, DOI `10.1016/j.cels.2022.05.002`.

The earlier v2.1 GSE comparison code is reusable only for authoritative H5AD loading, `layers['counts']`, gene mapping, raw-count aggregation, and frozen-parser application. Its binary treatment recoding, merging of the two control arms by numeric animal label, and paired-delta logic are not reusable in v3.1.

Cells are not independent validation samples.

### 6.2 Locked analyses

1. Coverage audit before scores are interpreted.
2. Baseline young-versus-aged direction, medians, median difference, and M1-M4 agreement are calculated separately in `Tg+/Dox-` and `Tg-/Dox+` controls. Direction must agree in both control arms before baseline direction is called preserved. A pooled-control result is secondary and retains each arm-specific animal as a separate unit.
3. Within each age, estimate two unpaired SOKM contrasts:

$$
\Delta_{Tg+/Dox-}=\operatorname{center}(Y_{Tg+/Dox+})-
\operatorname{center}(Y_{Tg+/Dox-}),
$$

and:

$$
\Delta_{Tg-/Dox+}=\operatorname{center}(Y_{Tg+/Dox+})-
\operatorname{center}(Y_{Tg-/Dox+}).
$$

4. Report all arm-specific animal scores, mean and median contrasts, and descriptive uncertainty. A SOKM direction is called concordant only when its sign agrees against both control arms.
5. Analyze young and aged strata separately before any pooled summary.
6. Audit cell-state composition by nested animal and exact treatment arm.
7. Perform state-restricted unpaired treatment-arm comparisons only for states represented in SOKM and the relevant control arm. Reprogramming states absent from controls are composition findings, not within-state contrasts.
8. Compare M1 and M4 treatment contrasts directly. Do not calculate a cross-treatment per-animal delta correlation because animals are not paired across treatment arms.
9. Relate Youth contrasts to prespecified MSC identity, ECM/collagen, proliferation/cell-cycle, and reprogramming programs only if their gene sets are frozen before Youth outcomes are inspected.

GSE176206 provides binary young/aged labels rather than a continuous age series. It therefore does not test continuous age ordering.

The median contrast is the primary descriptive treatment effect; the mean contrast is supporting. Both are reported, but model roles are not selected from either magnitude alone.

### 6.3 Small-sample rules

There are only three known animals per age-treatment arm. Therefore:

- all arm-specific animal values must be shown;
- there is no cross-treatment animal pairing;
- correlations within an arm based on three animals are descriptive;
- bootstrap intervals, if reported, are descriptive and cannot establish model superiority;
- no cell-level p-value is used as donor-level evidence;
- a state-restricted arm comparison requires at least two animals per represented arm and remains descriptive;
- missing state-condition combinations are reported as not assessable, not treated as zero.

### 6.4 Identity and Risk boundary

No formally frozen standalone Identity or Risk scorer has been identified in the current Youth Score workspace. Therefore formal `Youth-Identity` and `Youth-Risk` model coupling are currently `not_assessable`.

Existing candidate marker modules may be used only as clearly labeled exploratory program scores after their definitions are frozen. They must not be described as validated Identity or Risk models, and they cannot trigger Youth model changes.

### 6.5 Interpretation

- concordant baseline direction supports external directional consistency;
- concordant unpaired SOKM contrasts against both controls support shared perturbation sensitivity;
- neither result alone proves rejuvenation or accelerated aging;
- Youth decrease together with loss of MSC identity programs is more consistent with state or identity remodeling than with a simple aging trajectory;
- M1-M4 disagreement is reported as biological information and is not resolved by retuning.

## 7. Phase 2: TMS Droplet Cross-Assay Transportability

Droplet is run after GSE176206. It is a within-TMS cross-assay transportability and robustness check, not fully independent external validation and not a winner-selection dataset.

### 7.1 Input and units

Use the previously audited processed inputs:

- `data_droplet/processed/tms_limb_msc_pseudobulk_counts.rds`;
- `data_droplet/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv`.

The frozen input contract is 20,138 genes by 12 mice, with metadata in exactly the same mouse order. The age counts are two 3-month, four 18-month, two 21-month, and four 24-month mice. The processed objects were derived from the cleaned Limb Muscle MSC data and are the same source objects used by the v2.1 FACS-to-Droplet comparison.

Both models must receive this identical mouse-level raw-count pseudobulk matrix and metadata. No Droplet outcome may alter a frozen artifact.

For M1/M4 scoring, the raw all-gene denominator is recomputed as `colSums(tms_limb_msc_pseudobulk_counts.rds)` and must equal metadata `total_counts`. Metadata `raw_library_size` is a legacy 12,735-gene filtered-DGE library size and must be renamed `filtered_DGE_library_size` in analysis outputs; it must not be passed to the frozen parser as the all-gene denominator. Metadata `effective_library_size` and `data_droplet/processed/pseudobulk_dge_tmm.rds` are used only for technical-association auditing. `pseudobulk_logcpm.rds` is not an M1/M4 scoring input.

### 7.2 Locked analyses

For M1 and M4 report:

- coverage and missing genes;
- donor scores and paired M4-minus-M1 score differences;
- young and old medians and young-minus-old median difference;
- AUC and all-age Spearman correlation;
- old-only correlation and ordering across 18, 21, and 24 months as descriptive sensitivities;
- correlation with all-gene library size, filtered-DGE library size, TMM effective library size, cell count, and detected genes;
- overall, young-only, old-only, and age-group-residualized M1-M4 agreement;
- donor rank differences and the largest disagreements.

Use paired donor bootstrap on the same resampled mice for differences in AUC, all-age correlation, old-only correlation, median separation, and technical correlations. These intervals quantify uncertainty and are not tuning criteria.

### 7.3 Interpretation

- direction preserved by both models: cross-assay directional transportability;
- M1 preserved and M4 degraded: evidence that stricter stability may be FACS-specific;
- M4 preserves stronger ordering without worse coverage, technical association, or uncertainty: support for a stronger external role, not proof of universal superiority;
- both retain direction but remain technically associated: mixed transportability without technical independence;
- both fail direction with adequate coverage: failed transportability under this assay/cohort shift;
- inadequate coverage: feature-limited application, not evidence of biological failure.

## 8. Core M1-M4 Agreement Analysis

Overall correlation is insufficient because both models may merely separate young from old. In each dataset report, where sample size permits:

$$
\rho_{overall},\qquad
\rho_{young},\qquad
\rho_{old},\qquad
\rho_{age-adjusted}.
$$

For GSE176206, additionally report within-exact-treatment-arm agreement and agreement after removing age-by-exact-treatment centers. Compare M1 and M4 arm-level treatment contrasts directly. Do not report a per-animal treatment-delta correlation because animal labels are nested within treatment arms.

Interpretation is frozen:

- high overall and adjusted agreement: largely shared latent axis;
- high overall but low adjusted agreement: concordance mainly reflects group separation;
- low overall and adjusted agreement: potentially distinct biological axes.

## 9. Final Role Assignment

No single point estimate or p-value selects a winner. Final roles use the joint evidence from the completed permutation, GSE biological interpretation, Droplet transportability, coverage, uncertainty, technical association, and M1-M4 adjusted agreement.

### Scenario A: M1 and M4 are externally equivalent

Retain M1 as the main deployable reference; retain M4 as sensitivity comparator.

### Scenario B: M4 is consistently more informative externally

M4 may be upgraded to the preferred high-stringency candidate only if it adds interpretable GSE information and Droplet transportability without meaningful degradation of coverage, technical association, or uncertainty. Its post-ablation status remains explicit.

### Scenario C: M4 improves internally but degrades externally

Treat its internal trajectory improvement as cohort-specific optimization and retain M1 as the deployable reference.

### Scenario D: neither model transports or is interpretable

Report the negative result. Do not tune v3.1 on external outcomes. Any new model requires a separately frozen development protocol.

## 10. Stop Rules

After this freeze:

- do not tune genes, weights, thresholds, normalization, calibration, orientation, or coverage gates;
- do not remove inconvenient donors or animals;
- do not add M2, M3, or a new model after inspecting outcomes;
- do not use unknown GSE animals in primary analyses;
- do not pair equal numeric GSE animal labels across treatment arms;
- do not merge `Tg+/Dox-` and `Tg-/Dox+` animals by numeric animal label;
- do not treat cells as independent replicates;
- do not use external data to recalibrate score scale;
- do not choose a model from p-value or effect magnitude alone;
- do not call Droplet independent external validation;
- do not call a GSE treatment response proof of rejuvenation, accelerated aging, or a chronological aging clock;
- do not invent Identity or Risk scores when frozen artifacts are unavailable.

## 11. Execution Order and Deliverables

The locked order is:

1. freeze and hash this protocol;
2. apply frozen M1/M4 to GSE176206;
3. complete coverage, baseline, unpaired treatment-arm, state, and exploratory program analyses;
4. apply the same frozen M1/M4 to TMS Droplet;
5. complete transportability, technical-association, agreement, and paired-bootstrap analyses;
6. assign model roles using the frozen scenarios;
7. write a combined final report without additional FACS tuning.

Required output roots:

- `outputs/facs_v3_1/external_validation/gse176206/`;
- `outputs/facs_v3_1/external_validation/droplet/`;
- `outputs/facs_v3_1/external_validation/final_synthesis/`.

Each phase must include input provenance, coverage tables, unit-level scores, uncertainty or descriptive limitations, figures, session information, script hashes, and a current manifest.

## 12. Claims Boundary

The intended final hierarchy is:

- FACS nested permutation: evidence that the frozen training pipeline contains non-random age-label association;
- GSE176206: external directional and biological-perturbation interpretability;
- TMS Droplet: cross-assay directional transportability and robustness;
- combined evidence: assignment of model roles, not proof of a universal aging clock.
