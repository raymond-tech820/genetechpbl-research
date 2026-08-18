# A Donor-Level MSC Youth Score Framework Reveals Trade-offs in Small-Cohort Aging-State Modeling

## Abstract

Single-cell transcriptomic studies of mesenchymal stromal/stem cell (MSC) aging often contain many cells but few independent donors, creating a tension between apparent predictive performance and donor-level generalizability. We developed a donor-level Youth Score framework that aggregates raw single-cell counts into mouse pseudobulk profiles, estimates sex-compatible age effects with factorial differential-expression models, and uses nested stability selection to construct weighted young-high and old-high gene modules. Fully nested leave-one-mouse-out (LOMO) analysis was used so that filtering, normalization, feature selection, weighting, standardization, and calibration were rebuilt without the held-out mouse.

Controlled ablations showed that age discrimination, within-old age ordering, technical association, and fold-level signature stability reached their strongest values in different model branches. We therefore froze two scientifically distinct models: M1, a 58-gene primary operational reference selected at an inner-fold frequency threshold of $\pi_g\geq0.75$, and M4, a 29-gene post-ablation high-stringency sensitivity model selected at $\pi_g\geq0.90$. Both achieved nested LOMO AUC 1.000. Their all-age Spearman correlations were $-0.806$ and $-0.844$, respectively, while correlations with raw library size remained substantial ($-0.692$ and $-0.727$). In 999 fully nested, sex-stratified age-label permutations, the prespecified absolute all-age correlations exceeded the randomized-label null for M1 ($p=0.023$) and M4 ($p=0.010$).

Frozen application to GSE176206 preserved the young-higher-than-aged direction in two control arms and produced negative median SOKM-minus-control contrasts in both age groups. Frozen application to TMS Droplet preserved young-old separation, negative age ordering, and complete signature coverage across assay. M1 and M4 remained highly concordant in GSE176206 and Droplet, and paired Droplet bootstrap intervals did not establish a performance difference between them. These results support a rigorously characterized MSC youth-associated transcriptional-state framework. They also expose a central small-cohort constraint: stronger descriptive trajectory ordering can accompany stronger technical coupling and lower fold-level signature overlap.

## 1. Introduction

### 1.1 Biological and statistical problem

Quantifying a young-like MSC transcriptional state could support studies of aging, regenerative interventions, and partial reprogramming. The available single-cell datasets create a specific statistical challenge: thousands of cellular measurements may arise from only a small number of biologically independent animals. Treating cells as replicates inflates effective sample size and allows donor-specific, technical, or cohort-specific structure to masquerade as generalizable age biology.

Partial-reprogramming experiments further show why youth-associated expression and cell identity require separate interpretation. Roux et al. reported partial-reprogramming-associated restoration of youthful expression programs together with transient suppression of somatic identity programs in murine systems. This biological separation motivates a Youth-state score whose response is interpreted alongside independent identity and risk axes (Roux et al., 2022).

Initial exploratory work showed that three objectives require separate evaluation in these cohorts:

1. **young-old discrimination**, or separation of broad cohort states;
2. **age ordering**, including ordering among older donors;
3. **robustness**, covering technical association and feature-selection stability.

The final framework addresses these objectives at the donor level and treats technical association as an explicit model diagnostic. Earlier Droplet and FACS models supplied feasibility evidence and exposed donor imbalance, sex-age structure, assay depth, and subtissue contamination as consequential sources of uncertainty. The manuscript analysis begins from the cleaned FACS cohort and uses earlier model behavior only as motivation for prespecified ablations.

### 1.2 Contribution

We establish a donor-level MSC Youth Score framework with three linked contributions:

1. donor-level pseudobulk scoring with sex-compatible factorial age-effect estimation;
2. fully nested stability-controlled feature selection and explicit characterization of the trade-off among age association, technical burden, and signature stability;
3. frozen evaluation under regenerative perturbation and cross-assay transfer.

The central methodological finding is a Pareto trade-off under small-cohort constraints. Increasing stability-stringency improved descriptive age ordering in this dataset, while raw-library association strengthened and cross-fold signature Jaccard decreased. The final output is therefore a frozen scoring framework with a primary operational model and a high-stringency sensitivity model, each assigned a distinct evidential role.

## 2. Methods

### 2.1 Training cohort and metadata audit

The final training cohort contained 815 FACS-profiled TMS limb-muscle MSCs from 14 mice after removal of 120 diaphragm cells that were confined to old animals. All retained cells had standardized forelimb/hindlimb subtissue labels. The expression matrix contained 22,966 genes.

| Age | Female mice | Male mice | Female cells | Male cells |
|---|---:|---:|---:|---:|
| 3 months | 2 | 4 | 130 | 134 |
| 18 months | 2 | 2 | 158 | 76 |
| 24 months | 0 | 4 | 0 | 317 |

The 24-month group contains only male mice. The factorial design supports female and male young-old contrasts, while a female 18-to-24-month trajectory is unidentifiable. Four outer folds were retained with explicit weak-support flags because one factorial cell had only one training mouse. Matrix counts were used as the expression source of truth; two small metadata `n_counts` discrepancies in donor `24_59_M` did not alter the pseudobulk matrix or model output.

### 2.2 Donor-level pseudobulk

Raw single-cell counts were summed within each mouse:

$$
C_{gm}=\sum_{c\in m} C_{gc},
$$

where $C_{gc}$ is the raw count for gene $g$ in cell $c$, and $C_{gm}$ is the mouse-level pseudobulk count. Cells are measurement units; mice are inference and validation units.

This sum-based pseudobulk preserves the count scale required by edgeR TMM and limma-voom. It also prevents donor cell abundance from being converted into artificial biological replication.

### 2.3 Training-partition filtering, TMM, and voom

Within every training partition, gene $g$ was retained when

$$
\mathrm{CPM}_{gm}>1
$$

in at least two training mice. Outer LOMO folds retained 12,535 to 12,777 genes, with a median of 12,708. TMM factors and the voom mean-variance relationship were estimated from the training mice only. TMM and voom were used for differential-expression estimation; deployment scoring used the frozen raw all-gene transform described below.

PCA was used as exploratory quality control for age, sex, library size, detected genes, and cell count. It did not determine signature membership.

### 2.4 Sex-compatible age-effect estimation

An age-only model supplied a direction-consistency check:

$$
E_g\sim \mathrm{AgeGroup}.
$$

The main differential-expression model was a factorial cell-means design:

$$
E_g\sim 0+\mathrm{Sex}:\mathrm{AgeGroup}.
$$

For each gene, female, male, common, and interaction effects were defined as

$$
\beta_{g,F}=\mu_{g,F,O}-\mu_{g,F,Y},
$$

$$
\beta_{g,M}=\mu_{g,M,O}-\mu_{g,M,Y},
$$

$$
\beta_{g,C}=\frac{\beta_{g,F}+\beta_{g,M}}{2},
$$

$$
\beta_{g,I}=\beta_{g,M}-\beta_{g,F}.
$$

Genes entered the reliability pool when female and male effects shared a nonzero direction, the age-only and common factorial effects shared a nonzero direction, common statistics were finite, and the gene was absent from the frozen sex-linked exclusion list. This construction emphasizes age effects that remain directionally compatible across sex despite the unbalanced age-by-sex cohort.

### 2.5 Interaction-aware ranking

Sex heterogeneity was summarized by

$$
R_g=\frac{|\beta_{g,I}|}{|\beta_{g,F}|+|\beta_{g,M}|+\epsilon},
$$

with interaction penalty

$$
P_g=\frac{1}{1+R_g}.
$$

Eligible genes were ranked by

$$
w_g=|\beta_{g,C}|\,|t_{g,C}|\,P_g,
$$

where $t_{g,C}$ is the moderated common-effect statistic. Genes with $\beta_{g,C}<0$ formed the young-high module; genes with $\beta_{g,C}>0$ formed the old-high module.

### 2.6 Nested stability selection

Inside each outer-training partition, provisional signatures were reconstructed in inner LOMO folds, with at most 50 genes per direction in each provisional signature. This cap stabilized inner-fold ranking under the very small donor count and was fixed before external application; it was not tuned against external performance. Gene-level selection frequency was

$$
\pi_g=
\frac{\text{eligible inner folds selecting }g}
{\text{number of eligible inner folds}}.
$$

M1 retained genes with $\pi_g\geq0.75$. M4 retained genes with $\pi_g\geq0.90$. The stable pool was reranked on the complete outer-training partition, capped at 100 genes per module, and never backfilled with genes below the frozen frequency threshold. Rank-deficient inner splits were recorded and skipped without changing the factorial model.

### 2.7 Frozen scoring function

M1 and M4 use the raw all-gene library denominator:

$$
X_{gm}=
\log_2\left(
10^6\frac{C_{gm}}{\sum_{h\in G_{all}}C_{hm}}+1
\right).
$$

Gene values were standardized using training-only parameters:

$$
Z_{gm}=\frac{X_{gm}-\mu_g^{train}}{s_g^{train}}.
$$

For module $A$, the weighted module score was

$$
S_A(m)=
\frac{\sum_{g\in A}w_gZ_{gm}}
{\sum_{g\in A}w_g}.
$$

The raw Youth Score was

$$
Y_m^{raw}=S_{young}(m)-S_{old}(m).
$$

Let $M_Y$ and $M_O$ be the median raw scores among young and old training mice. The calibrated score was

$$
Y_m=
\frac{Y_m^{raw}-M_O}{M_Y-M_O}.
$$

Training medians are fixed to 1 for young and 0 for old by construction. Model evaluation therefore uses held-out LOMO scores or frozen external applications.

### 2.8 Fully nested donor-level validation

For every held-out mouse, the complete training-side pipeline was rebuilt:

1. gene filtering;
2. TMM and voom estimation;
3. age-only and factorial differential expression;
4. reliability screening and interaction-aware ranking;
5. inner LOMO stability selection;
6. module selection and weight estimation;
7. training-only standardization and calibration;
8. frozen-transform scoring of the untouched held-out mouse.

This design prevents feature-selection, normalization, and calibration leakage. Performance was summarized by young-old AUC, all-age Spearman correlation, old-only Spearman correlation, young-old median difference, associations with technical variables, signature size, and fold Jaccard.

### 2.9 Controlled ablation design and operating points

The exact historical factorial stability-selected implementation was retained as M0, a sanity anchor for testing whether later changes altered the established donor ordering. Three independent design questions were then evaluated:

| Analysis | Role | Design question |
|---|---|---|
| M1 | Primary operating point | Does correcting the deployment denominator preserve the historical solution? |
| M2 | Technical-adjustment ablation | Does adding library size and cell count during DE reduce final score-level technical coupling? |
| M3 | Stability-threshold sweep | How does $\pi_g=0.50$, 0.75, or 0.90 move the score across age ordering, technical association, and fold stability? |
| M4 | High-stringency operating point | What behavior is retained when the corrected denominator is combined with $\pi_g\geq0.90$? |

M1, M2, and M3 branch independently from M0. M2 retains the M0 normalization contract, so its effect isolates the DE adjustment. M4 was frozen after review of the independent ablations and retains a post-ablation exploratory role. Complete branch definitions and results are reported in Appendices A-C.

### 2.10 Frozen permutation test

The primary permutation statistic was

$$
T_{\rho}=|\rho(Y,\mathrm{Age})|.
$$

Complete 3-, 18-, and 24-month labels were reassigned within sex strata. Every assignment reran the fully nested pipeline, including feature selection and calibration. The empirical two-sided p-value was

$$
p=\frac{1+\sum_{b=1}^{999}I(T_b\geq T_{obs})}{1000}.
$$

All 999 assignments were unique and valid for both models. Young-old AUC distance from 0.5 was prespecified as supporting evidence; young-old median difference was recorded as a non-primary statistic.

### 2.11 Frozen external application

External datasets were never used to change genes, weights, standardization, calibration, thresholds, or score orientation.

For GSE176206, the MSC SOKM object contained 20,661 source cells, of which 18,686 met the frozen analysis criteria and formed 18 age-treatment-animal units. Raw counts from `layers['counts']` were aggregated by age, exact treatment arm, and nested animal identifier. Numeric animal labels were treated as nested within each exact arm, so treatment comparisons were unpaired. The source study profiled young and aged murine MSCs after transient SOKM induction and chase; the present analysis used this dataset to assess baseline age direction and perturbation-related biological consistency (Roux et al., 2022; National Center for Biotechnology Information, 2022).

For TMS Droplet, the frozen models were applied to 12 mouse-level pseudobulk profiles: two at 3 months, four at 18 months, two at 21 months, and four at 24 months. The 20,138-gene all-gene column sum supplied the parser denominator. This analysis assessed within-TMS cross-assay transportability. Paired donor bootstrap compared M4-minus-M1 differences without model reselection.

## 3. Results

### 3.1 Small-cohort evaluation exposes multiple optimization objectives

Controlled ablations were designed to identify representative operating points across competing objectives. The diagnostic branches separated broad age discrimination from within-old ordering, technical association, and signature stability. Every model retained 14 valid held-out scores and zero empty outer signatures.

M0 served as a provenance anchor, establishing the reproducible v2.1 solution as the fixed baseline for every subsequent v3.1 modification.

| Analysis | Main change | Key observation |
|---|---|---|
| M1 | Corrected raw all-gene denominator at $\pi_g\geq0.75$ | Preserved the historical donor ranking and supplied the primary operational reference |
| M2 | Added library size and cell count during DE | Produced limited raw-library improvement while weakening age ordering and increasing cell-count coupling |
| M3 | Swept $\pi_g$ from 0.50 to 0.90 | Showed that higher stringency strengthened descriptive age ordering while increasing raw-library coupling and lowering fold Jaccard |
| M4 | Combined the corrected denominator with $\pi_g\geq0.90$ | Defined a representative high-stringency operating point for sensitivity analysis |

M1 preserved the M0 donor ranks after correcting the deployment denominator. This correction separated scoring normalization from training-fold gene filtering without changing the learned donor ordering. M2 supplied the central technical-adjustment control: raw-library correlation improved only modestly, while old-only ordering changed from $-0.109$ to $+0.436$ and cell-count association strengthened. The M3 threshold sweep then identified stability stringency as the design factor that moved the score toward stronger descriptive age ordering, accompanied by stronger raw-library coupling and lower fold stability.

![Stability-threshold trade-off under nested LOMO](outputs/facs_v3_1/controlled_ablation/figures/manuscript_stability_tradeoff.png)

Together, these ablations explain the selection of M1 and M4 as representative operating points. Complete M0-M4 metrics, denominator-equivalence results, and threshold-sweep values are reported in Appendices A-C.

### 3.2 Two frozen models define representative operating points

Full-data freezing selected 58 genes for M1, comprising 30 young-high and 28 old-high genes, and 29 genes for M4, comprising 16 young-high and 13 old-high genes. Calibration denominators were 2.65029 and 2.70224.

| Model | Frozen role | Stability threshold | Genes | Nested AUC | All-age $\rho$ | Old-only $\rho$ |
|---|---|---:|---:|---:|---:|---:|
| M1 | Primary operational reference | 0.75 | 58 | 1.000 | -0.806 | -0.109 |
| M4 | High-stringency sensitivity model | 0.90 | 29 | 1.000 | -0.844 | -0.327 |

M1 preserves the corrected deployment transform and the broader stable pool. M4 emphasizes genes selected in at least 90% of eligible inner folds and provides a trajectory-sensitive operating point. M4's development followed review of the ablations, so its favorable internal point estimates remain exploratory evidence.

M4 was retained to quantify how conclusions change under a stricter stability constraint. This sensitivity role tests robustness to feature-stringency while preserving M1 as the primary operational reference.

![Nested LOMO scores for the frozen candidates](deliverables/facs_msc_youth_score_v3_1/training_experiment/training_nested_lomo_scores_by_age.png)

![Frozen signature composition](deliverables/facs_msc_youth_score_v3_1/training_experiment/frozen_signature_module_composition.png)

### 3.3 Frozen nested pipelines exceed randomized-age nulls

The 999 sex-stratified permutations reconstructed every model component within each outer fold. All assignments produced 14 valid scores per model, preserved sex-age counts, and had no zero-signature or rank-deficient outer folds.

| Model | Statistic | Observed | Null mean | Null SD | $Z_{null}$ | Empirical $p$ |
|---|---|---:|---:|---:|---:|---:|
| M1 | $\rho_{age}$ | 0.806 | 0.358 | 0.229 | 1.953 | 0.023 |
| M1 | $\mathrm{AUC}-0.5$ | 0.500 | 0.218 | 0.138 | 2.036 | 0.039 |
| M4 | $\rho_{age}$ | 0.844 | 0.355 | 0.225 | 2.177 | 0.010 |
| M4 | $\mathrm{AUC}-0.5$ | 0.500 | 0.215 | 0.139 | 2.047 | 0.036 |

The prespecified result supports age-label association beyond the practical randomized-label null for both frozen pipelines. The test addresses label association under the cohort structure; technical independence and continuous-clock validity require separate evidence.

![Formal permutation null for absolute all-age correlation](deliverables/facs_msc_youth_score_v3_1/training_experiment/formal_permutation_abs_all_age_rho_null.png)

### 3.4 Frozen evaluation supports biological direction and assay transportability

#### 3.4.1 GSE176206 biological application

Both models passed total and module-level coverage gates. Minimum module gene and weighted coverage were 0.867 and 0.870 for M1, and 0.938 and 0.922 for M4. In both control arms, median Youth Scores were higher in young than aged animals.

| Model | Control arm | Young median | Aged median | Young-minus-aged median | AUC |
|---|---|---:|---:|---:|---:|
| M1 | `Tg+/Dox-` | 0.4864 | 0.4298 | 0.0566 | 1.000 |
| M1 | `Tg-/Dox+` | 0.5003 | 0.4550 | 0.0453 | 0.778 |
| M4 | `Tg+/Dox-` | 0.4137 | 0.3635 | 0.0503 | 1.000 |
| M4 | `Tg-/Dox+` | 0.4481 | 0.3856 | 0.0626 | 0.778 |

The prespecified median SOKM-minus-control contrast was negative for both models, both ages, and both control arms.

| Model | Age | Versus `Tg+/Dox-` | Versus `Tg-/Dox+` |
|---|---|---:|---:|
| M1 | Young | -0.0648 | -0.0787 |
| M1 | Aged | -0.0148 | -0.0400 |
| M4 | Young | -0.0806 | -0.1150 |
| M4 | Aged | -0.0455 | -0.0676 |

The M1 aged SOKM-versus-`Tg+/Dox-` mean difference was slightly positive at 0.0042 while the prespecified median difference was -0.0148. The perturbation result is therefore directional median evidence across small unpaired arms, with heterogeneous support across summary statistics.

M1 and M4 had overall Spearman agreement 0.930 and age-treatment-residualized agreement 0.907 across the 18 nested-animal treatment units. State-restricted contrasts were heterogeneous, and pure reprogramming states absent from controls represented composition changes. These observations localize the strongest evidence to an aggregate perturbation-sensitive youth-associated state axis. Formal separation of Youth, MSC identity, and treatment risk awaits independently frozen Identity and Risk scorers.

![GSE176206 scores by exact treatment](deliverables/facs_msc_youth_score_v3_1/cross_assay_comparison/gse176206_scores_by_exact_treatment.png)

#### 3.4.2 TMS Droplet cross-assay transportability

All M1 and M4 signature genes were present in the Droplet pseudobulk matrix. Both frozen models achieved AUC 1.000, young medians above old medians, and negative all-age and old-only correlations.

| Model | Young median | Old median | AUC | All-age $\rho$ | Old-only $\rho$ | All-gene library $\rho$ | Cell-count $\rho$ | Detected-gene $\rho$ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| M1 | 0.4455 | 0.3152 | 1.000 | -0.808 | -0.662 | -0.028 | -0.098 | 0.336 |
| M4 | 0.3882 | 0.2547 | 1.000 | -0.830 | -0.701 | -0.035 | -0.105 | 0.308 |

M1 and M4 donor scores had overall Spearman agreement 0.965, old-only agreement 0.939, and exact-age-residualized agreement 0.846. M4 had slightly more negative age-correlation point estimates, while every prespecified paired-bootstrap interval for M4-minus-M1 included zero. The Droplet analysis therefore supports shared cross-assay direction and leaves the frozen model roles unchanged.

![Droplet scores by age](deliverables/facs_msc_youth_score_v3_1/cross_assay_comparison/droplet_scores_by_age.png)

### 3.5 Frozen parser reproduces the training implementation

The exported parsers reproduced independently computed training scores with maximum calibrated differences of $2.00\times10^{-15}$ for M1 and $1.18\times10^{-15}$ for M4. Adding a new sample left existing sample scores unchanged, full-coverage inputs passed, and a deliberately feature-limited input triggered the frozen coverage error.

These controls establish numerical reproducibility and deployable normalization. They complement the inferential analyses without using apparent full-data scores as performance evidence.

## 4. Discussion

### 4.1 A youth-associated state score under small-cohort constraints

The framework derives a reproducible donor-level youth-associated transcriptional axis within the evaluated MSC cohorts. Its key advance is the integration of sex-compatible effect estimation, inner-fold feature stability, fully nested donor validation, technical diagnostics, and frozen external application. This design addresses the gap between large cell counts and small independent-donor counts that characterizes many single-cell aging datasets.

The formal permutation results show that the frozen pipelines contain stronger age-label association than expected under sex-stratified randomized labels. GSE176206 and Droplet extend this evidence by showing directionally consistent behavior under a regenerative perturbation context and an assay shift. Together, these layers support a characterized Youth-state framework through complementary statistical, biological, and transportability evidence.

### 4.2 Biological sensitivity, technical robustness, and stability form a trade-off

The controlled ablations provide the main methodological result. Technical-covariate adjustment during DE changed selected genes substantially and did not deliver coherent score-level technical improvement. Lower stability stringency reduced some technical correlations while weakening age behavior. Higher stringency strengthened descriptive age ordering while increasing raw-library coupling and reducing fold Jaccard.

This finding has two implications. First, technical correction at one stage of a pipeline does not guarantee technical independence of the final score. Second, stricter gene-selection frequency can enrich a smaller, more age-ordered signature while reducing cross-fold gene overlap because few genes satisfy the threshold consistently across all small training partitions. The appropriate operating point therefore depends on the biological and deployment question.

### 4.3 Model roles follow the evidence hierarchy

M1 remains the primary operational reference because it was frozen as the corrected-denominator reference before external application and retains a broader stable signature. M4 remains a high-stringency sensitivity model because it was selected after inspecting the ablation results. The two models are highly concordant in both external applications, and Droplet paired-bootstrap comparisons provide no evidence of a material performance difference.

This role-based interpretation preserves the information in M4 without allowing post-ablation internal gains to overwrite selection history. It also gives future studies a direct sensitivity analysis: conclusions that survive both M1 and M4 are less dependent on the stability threshold.

### 4.4 Youth Score interpretation in regenerative experiments

In GSE176206, SOKM exposure shifted median scores downward in young and aged groups. The same perturbation created heterogeneous state-restricted effects and states absent from controls. The Youth Score therefore responds to perturbation-associated transcriptional remodeling, with aggregate changes potentially reflecting both within-state expression shifts and state-composition changes.

Roux et al. originally reported restoration of their youthful-expression measures together with transient loss of somatic identity. The downward shift of the frozen limb-MSC Youth Score highlights that youth-associated scores can emphasize different transcriptional axes. This contrast strengthens the need to interpret the present score as a cohort-trained state axis and to evaluate it jointly with independently defined Identity and Risk measures (Roux et al., 2022).

A favorable regenerative interpretation requires a joint state across complementary axes:

$$
Y_{Youth}\uparrow
\land I_{Identity}\ \mathrm{preserved}
\land R_{Risk}\ \mathrm{low}.
$$

Independent Identity and Risk models were unavailable before the frozen GSE analysis. Their coupling with Youth remains a defined future analysis, and no ad hoc marker proxy was introduced.

### 4.5 Limitations

1. The FACS training cohort contains 14 mice, including only two female and two male mice at 18 months and no females at 24 months.
2. Four outer folds contain weak factorial support, although all outer designs remained full rank and all held-out scores were finite.
3. FACS nested scores remain correlated with raw library size and cell count; technical independence has not been demonstrated.
4. M4 was selected after reviewing stability-threshold ablations and therefore retains an exploratory role even after a favorable permutation result.
5. GSE176206 contains three animals per age-treatment arm and treatment arms are unpaired; effect estimates and bootstrap intervals are descriptive.
6. TMS Droplet provides a within-study cross-assay assessment, with biological and annotation relationships to the FACS cohort.
7. The available datasets support young-old state separation more strongly than a continuous chronological aging trajectory.

### 4.6 Conclusion

The Youth Score study establishes a donor-level MSC youth-associated scoring framework for small single-cell cohorts. Fully nested inference, factorial age modeling, stability-controlled feature selection, and frozen external application produced a reproducible young-like state axis. Controlled ablations showed that descriptive trajectory sensitivity, technical robustness, and signature stability occupy different operating points within the tested design space. M1 provides the primary operational score, and M4 supplies a frozen post-ablation high-stringency sensitivity analysis for downstream studies.

## 5. Appendices

### Appendix A. Complete M0-M4 ablation results

All values below use held-out mouse scores from fully nested LOMO. They are retained as supporting ablation evidence and are not apparent full-data fits.

| Branch | AUC | All-age $\rho$ | Old-only $\rho$ | Young-old median difference | Raw-library $\rho$ | Effective-library $\rho$ | Cell-count $\rho$ | Detected-gene $\rho$ | Median genes | Median fold Jaccard |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| M0 exact reference | 1.000 | -0.806 | -0.109 | 0.2602 | -0.692 | -0.473 | -0.557 | +0.042 | 54.5 | 0.407 |
| M1 raw denominator | 1.000 | -0.806 | -0.109 | 0.2601 | -0.692 | -0.473 | -0.557 | +0.042 | 54.5 | 0.407 |
| M2 adjusted DE | 1.000 | -0.713 | +0.436 | 0.2677 | -0.653 | -0.508 | -0.590 | +0.037 | 52.0 | 0.422 |
| M3, $\pi=0.50$ | 0.979 | -0.670 | +0.436 | 0.3029 | -0.552 | -0.459 | -0.506 | -0.095 | 77.0 | 0.417 |
| M3, $\pi=0.75$ | 1.000 | -0.806 | -0.109 | 0.2602 | -0.692 | -0.473 | -0.557 | +0.042 | 54.5 | 0.407 |
| M3, $\pi=0.90$ | 1.000 | -0.844 | -0.327 | 0.2729 | -0.727 | -0.468 | -0.559 | +0.099 | 31.0 | 0.352 |
| M4 raw denominator, $\pi=0.90$ | 1.000 | -0.844 | -0.327 | 0.2728 | -0.727 | -0.468 | -0.559 | +0.099 | 31.0 | 0.352 |

Every branch had 14 valid outer scores, zero empty or invalid signatures, complete held-out gene coverage, and four weak-support folds.

### Appendix B. Historical reference and denominator equivalence

M0 exactly reproduced the historical v2.1 `factorial_stability_selected` pseudobulk, metadata, nested scores, fold signatures, diagnostics, and summary. M1 then changed only the deployment scoring denominator from the fold-filtered gene library to the raw all-gene library.

| Comparison | Spearman correlation | Median absolute score change | Maximum absolute score change | Signature identity |
|---|---:|---:|---:|---|
| M1 versus M0 | 1.000 | 0.000041 | 0.000158 | Exact |

The corrected denominator preserved every donor rank and every selected gene. This establishes M1 as a deployment-contract correction of the historical solution, with unchanged biological ordering.

### Appendix C. Stability-threshold sensitivity

M3 varied only the inner selection-frequency threshold while retaining the M0 DE design and historical scorer. The $\pi=0.75$ branch was an implementation control and reproduced M0 exactly.

| Threshold | AUC | All-age $\rho$ | Old-only $\rho$ | Raw-library $\rho$ | Cell-count $\rho$ | Median genes | Median fold Jaccard | Pooled gene Jaccard versus M0 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.50 | 0.979 | -0.670 | +0.436 | -0.552 | -0.506 | 77.0 | 0.417 | 0.755 |
| 0.75 | 1.000 | -0.806 | -0.109 | -0.692 | -0.557 | 54.5 | 0.407 | 1.000 |
| 0.90 | 1.000 | -0.844 | -0.327 | -0.727 | -0.559 | 31.0 | 0.352 | 0.689 |

M4 transferred the $\pi=0.90$ signature to the corrected raw all-gene denominator. M4 and M3 $\pi=0.90$ had score Spearman correlation 1.000, median absolute calibrated difference $4.36\times10^{-5}$, and maximum difference $1.74\times10^{-4}$. The threshold change generated the shift toward stronger descriptive age ordering; the denominator correction preserved that ordering.

The trade-off figure in Section 3.1 was generated directly from the locked branch summary.

### Appendix D. Additional prespecified permutation statistic

Young-old median difference was recorded before the formal run but was not designated as a primary or supporting test. It is retained here to preserve the complete permutation result without elevating its inferential role.

| Model | Observed young-old median difference | Null mean | Null SD | $Z_{null}$ | Empirical $p$ |
|---|---:|---:|---:|---:|---:|
| M1 | 0.2601 | -0.0767 | 0.1294 | 2.603 | 0.002 |
| M4 | 0.2728 | -0.0868 | 0.1621 | 2.219 | 0.007 |

### Appendix E. Model-development provenance

Earlier model stages are retained as methodological provenance:

- **Droplet v1** established feasibility and revealed the risks of sparse donor replication and cohort imbalance.
- **FACS v2/v2.1** introduced donor-richer pseudobulk modeling, sex-compatible factorial analysis, nested LOMO, low-depth and sex-sensitivity checks, and the reproducible stability-selected reference.
- **FACS v3** tested technical-aware ranking and showed that reduced technical association could accompany weaker age behavior.
- **FACS v3.1** converted these observations into controlled ablations, froze M1 and M4, and completed formal permutation and external application.

Detailed version-specific analyses, early parameter sweeps, invalidated pre-fix checkpoints, and implementation audits remain outside the manuscript's principal evidence hierarchy.

## 6. References

Roux, A. E., Zhang, C., Paw, J., Zavala-Solorio, J., Malahias, E., Vijay, T., Kolumam, G., Kenyon, C., & Kimmel, J. C. (2022). Diverse partial reprogramming strategies restore youthful gene expression and transiently suppress cell identity. *Cell Systems, 13*(7), 574-587.e11. [https://doi.org/10.1016/j.cels.2022.05.002](https://doi.org/10.1016/j.cels.2022.05.002)

National Center for Biotechnology Information. (2022). *GSE176206: Partial reprogramming restores youthful gene expression through transient suppression of cell identity* [Gene expression dataset]. Gene Expression Omnibus. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176206)
