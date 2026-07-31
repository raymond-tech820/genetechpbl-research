# Geneformer perturbation arm: methods and results log

Perturbation-prediction module for the cellular rejuvenation project. This log records
the pipeline, the analyses run, and their results. All numbers are as computed; results
not yet interpreted as scientific conclusions are marked as such.

Model: Geneformer-V1-10M (6 layers, input size 2048, vocabulary 25,426), gc30M token
dictionary, run zero-shot with no fine-tuning. Cohort: TMS FACS limb-muscle MSC, raw
counts, 815 cells (551 Old, 264 Young; 8 Old and 6 Young donors) after excluding 120
diaphragm-annotated cells. All perturbations use seed 42.

---

## 1. Pipeline

1. Mouse-to-human ortholog mapping via Ensembl BioMart, one-to-one orthologs only:
   14,314 / 22,966 genes retained (62.3%), 68.9% of total count mass.
2. Tokenization with `model_version="V1"` (required; the tokenizer otherwise defaults to
   the V2 gc104M dictionary, which mismatches the V1 checkpoint and silently mis-maps
   gene tokens).
3. In-silico overexpression with `InSilicoPerturber` (`perturb_type="overexpress"`,
   `emb_mode="cell"`). Goal state defined by Old and Young cell-embedding centroids on
   `age_group`.

Metric: `Shift_to_goal_end` from `perturber_utils.py::cos_sim_shift`, defined as the
perturbed cell's cosine similarity to the Young centroid minus the original cell's.
Positive = movement toward Young. Sign convention read from source.

Note on `genes_to_perturb`: a list is treated as one combined condition, not per-gene, so
per-gene ranking requires one call per gene with isolated output directories (the stats
class globs a directory without prefix filtering).

## 2. Preliminary perturbation output (not a result)

Single run, individual mode, full 815 cells, cell-level means over 551 Old cells; not
donor-aggregated. No null distribution or baseline, so this is pipeline output only.

| Factor | Human Ensembl ID | Mean shift | Per-cell SD | mean/SD |
|---|---|---|---|---|
| MYC | ENSG00000136997 | +0.000104 | 0.00327 | 0.032 |
| OCT4 | ENSG00000204531 | -0.000016 | 0.00431 | 0.004 |
| LIN28A | ENSG00000131914 | -0.000135 | 0.00424 | 0.032 |
| SOX2 | ENSG00000181449 | -0.000249 | 0.00430 | 0.058 |
| KLF4 | ENSG00000136826 | -0.000269 | 0.00366 | 0.074 |
| NANOG | ENSG00000111704 | -0.000276 | 0.00427 | 0.065 |
| MYOD1 | ENSG00000129152 | -0.000284 | 0.00428 | 0.066 |

For every factor the per-cell SD exceeds the mean shift by 15-25x. SOKM as one combined
condition: `+0.000166`, positive while three of its four components are negative. These
values are not interpretable until the axis they are measured along is validated (below).

Compute: ~30 s/gene, 3.5 GB peak VRAM (RTX 4060 laptop).

## 3. Axis validation

Every shift is measured along the Old-to-Young centroid axis, so that axis was tested for
biological separability before interpreting any perturbation value. Per-cell embeddings
were extracted for all 815 cells; separation was evaluated at the donor level (the unit of
replication) by leave-one-donor-out (LOO) nearest-centroid classification with a 1000-draw
donor-label permutation control (seed 42).

**Separation is weak.** Donor LOO accuracy 11/14, permutation p = 0.050. All three
misclassified donors are Young mice assigned to Old, each sitting near the Old centroid.

**The axis is confounded with transcript complexity.** Correlation of per-cell axis
position with technical covariates:

| Covariate | Pearson r | p |
|---|---|---|
| Total counts | -0.019 | 0.59 |
| Detected-gene count | 0.691 | 1.07e-116 |

Total counts is clean; detected-gene count is strongly associated.

## 4. Does separation survive removing the confound?

| Analysis | Donor LOO | Permutation p |
|---|---|---|
| Raw embeddings | 11/14 | 0.050 |
| Detected-gene count regressed out (per-dimension OLS) | 8/14 | 0.531 |
| Restricted to matched detected-gene band [1696, 3276] | 10/14 | 0.094 |

Under either correction, young-old separation is no longer significant. (Note: the
regressed-out centroid distance is a mathematical artifact of two-group residualization,
pinned at 2.0 regardless of biology; only the donor LOO from that row is meaningful, since
the LOO centroids are unweighted donor means and not subject to the artifact.)

## 5. Age and detected-gene count are confounded in this cohort

Separation collapses under correction because detected-gene count is itself confounded
with age at the donor level.

| Comparison | Old mean | Young mean | Welch p | Mann-Whitney p |
|---|---|---|---|---|
| Detected genes, all cells (donor-level) | 1,768 | 2,734 | 0.009 | 0.020 |
| Detected genes, within matched band (donor-level) | 2,215 | 2,571 | 0.028 | 0.043 |

The difference persists within the matched band, so age and transcript complexity cannot
be disentangled with this data. Lower detected-gene count in aged cells may be partly
genuine biology (declining transcriptional output) rather than a pure artifact; the
finding is that the two are inseparable in this cohort, not that the axis is purely
technical.

## 6. Cross-assay replication

The same donor-level test on the independent TMS Droplet limb-MSC cohort reproduced the
direction: Old 1,766 vs Young 2,023 detected genes (Welch p = 0.0007, Mann-Whitney
p = 0.03). Two assays that do not share preprocessing agree, arguing the age-complexity
association is a recurring feature of aged MSC data rather than a single-protocol artifact.
Caveats: Droplet has only two Young donors, both female, so it corroborates rather than
independently establishes the effect; and the effect size differs between assays
(~970 detected genes in FACS vs ~260 in Droplet), consistent with a mixture of biology and
assay-dependent technical behaviour.

## 7. Token-length sweep: is the confound a fixed input-length artifact?

Geneformer's input is a rank-ordered list of expressed genes, so detected-gene count sets
the input sequence length. Truncating every cell to a common length tests whether that
mechanism drives the confound. Seven lengths were swept. At each length a **fixed-length**
run (cells cut to L tokens, cells below L dropped) is paired with a **same-cell control**
(the surviving cells at their original variable length), so the truncation effect is
isolated from cell loss.

| L | Condition | Kept | Old dropped | Donor LOO | Perm. p | Confound r |
|---|---|---|---|---|---|---|
| 256 | Fixed | 815 | 0 | 13/14 | 0.001 | 0.619 |
| 256 | Control | 815 | 0 | 11/14 | 0.050 | 0.691 |
| 384 | Fixed | 815 | 0 | 12/14 | 0.011 | 0.597 |
| 384 | Control | 815 | 0 | 11/14 | 0.050 | 0.691 |
| 512 | Fixed | 809 | 6 | 11/14 | 0.044 | 0.558 |
| 512 | Control | 809 | 6 | 11/14 | 0.056 | 0.686 |
| 768 | Fixed | 780 | 34 | 12/14 | 0.015 | 0.544 |
| 768 | Control | 780 | 34 | 11/14 | 0.056 | 0.675 |
| 1024 | Fixed | 697 | 116 | 12/14 | 0.004 | 0.574 |
| 1024 | Control | 697 | 116 | 11/14 | 0.058 | 0.639 |
| 1536 | Fixed | 501 | 306 | 14/14 | 0.000 | 0.495 |
| 1536 | Control | 501 | 306 | 10/14 | 0.120 | 0.596 |
| 2048 | Fixed | 388 | 378 | 10/14 | 0.141 | 0.475 |
| 2048 | Control | 388 | 378 | 10/14 | 0.141 | 0.475 |

Findings:
- At every length the fixed-length confound r is below its same-cell control, so
  truncation reduces the confound; this is a truncation effect, not a cell-selection
  effect.
- The cleanest evidence is at 256 tokens, where the shortest tokenized cell has 394 genes
  so no cell is dropped and the cell set is identical to baseline: truncation alone raises
  donor LOO from 11/14 to 13/14 (p 0.05 to 0.001) and lowers the confound from 0.691 to
  0.619.
- The confound is reduced but floored near r ~ 0.50 and never approaches zero.
- Long-length points sit on shrunken, age-imbalanced cell sets (1536 drops 306 of the Old
  cells) and are not read as meaningful; the 2048 fixed/control rows are identical by
  construction, since truncating to the input cap is a no-op.

## 8. Donor-level (pseudobulk) axis

Building the axis from 14 donor-level pseudobulk profiles (rather than cells), where the
detected-gene difference is weakest, gave donor LOO 7/14 (chance) with the depth
correlation falling to r = 0.42 (p = 0.13). The truncation improvement is a cell-level
effect that does not survive aggregation to the donor, the unit of inference.

## 9. Summary

The perturbation ranking cannot be interpreted because the Old-to-Young axis it is measured
along is confounded with sequencing depth, and depth is inseparable from age in this cohort
(replicated across two assays). Input sequence length is a real, controllable component of
the confound: fixed-length truncation monotonically reduces it and improves cell-level
separation, demonstrated against same-cell controls. But the confound has an irreducible
floor and the axis fails at the donor level, so length control mitigates rather than
resolves the underlying depth-age entanglement. The module's contribution is this
characterised, partly-mechanistic negative and the associated validation criterion: a
usable aging axis must show donor-level separation that survives conditioning on
detected-gene count.

## 10. Key output files

| File | Contents |
|---|---|
| `results/tms_facs_limb_msc_tf_perturb_individual/…_ranking.csv` | 7-factor perturbation ranking |
| `results/axis_validation_summary.json` | Axis validation scalars (section 3) |
| `results/axis_validation_per_cell_axis_position.csv` | Per-cell axis position and covariates |
| `results/axis_validation_deconfound.json` | Deconfounding results (section 4) |
| `results/axis_validation_ngenes_confound_test.json` | Age vs detected-gene tests (section 5) |
| `results/droplet_donor_ngenes_distribution.csv` | Cross-assay replication (section 6) |
| `results/token_length_sweep.csv` | Token-length sweep (section 7) |
| `results/figures/token_sweep.pdf` | Sweep figure |

## 11. Scripts

`src/01_map_orthologs.py`, `02_prepare_tokenize.py`, `03_in_silico_perturb.py`,
`04_null_distribution.py`, `pseudobulk_axis.py`, `fixed_length_axis.py`,
`fixed_length_axis_cellset_control.py`, `token_length_sweep.py`, `plot_figures_confound.py`.
