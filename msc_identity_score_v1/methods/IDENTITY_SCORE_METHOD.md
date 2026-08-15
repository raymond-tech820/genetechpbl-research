# Identity Score v1: Frozen Method

## 1. Scope

Identity Score v1 quantifies expression of a literature-supported mouse MSC
identity program in individual cells. It is a knowledge-driven module score,
not a trained classifier. The score is designed as an independent biological
axis alongside Youth and Risk.

The frozen identity core contains 11 positively directed genes:

`Cd34`, `Cd44`, `Cxcl12`, `Eng`, `Itgav`, `Lepr`, `Ly6a`, `Nt5e`, `Pdgfra`,
`Prrx1`, and `Thy1`.

`Vcam1` was reviewed but excluded from the v1 core because the available
evidence was not sufficiently specific to mouse limb MSCs. Complete sources,
URLs, rationales, and exclusion reasons are in `identity_gene_set_v1.csv`.

The broader `ecm_collagen_program` is contextual information and is **not** a
component of Identity Score v1.

## 2. Input and preprocessing

The scorer expects raw, non-negative cell-by-gene counts. For the frozen
GSE176206 run, counts were read from `layers/counts`, and mouse symbols were
read from `var/gene_name`.

For each cell c and gene g, raw count y is transformed as:

```text
x_cg = log(1 + 10000 * y_cg / sum_h(y_ch))
```

No batch correction, z-scoring across the full dataset, or treatment-dependent
feature selection is applied.

## 3. Primary score

The primary score is a signed module mean corrected by expression-matched
background genes. Because all v1 identity genes have direction `+1`, the
identity component is the mean normalized expression of the 11 included genes.

Background construction is deterministic:

1. Mean gene expression is estimated from 3,694 known-animal control cells.
2. Genes are divided into 24 expression bins.
3. For every identity gene, one eligible background gene is sampled from the
   same bin.
4. All frozen signature genes and predefined technical-noise genes are excluded
   from the background universe.
5. This procedure is repeated for 100 control sets with seed `20260729`.

For cell c:

```text
Identity_primary(c)
  = mean expression of included identity genes
  - mean score across 100 matched background sets
```

Higher values indicate stronger identity-program expression relative to genes
with similar average abundance.

## 4. Rank-based sensitivity score

The sensitivity method ranks the expressed genes within each cell. Each
identity gene contributes its within-cell expression percentile, while an
unexpressed identity gene contributes zero. The score is divided by the 11
frozen genes:

```text
Identity_rank(c)
  = sum of within-cell expression percentiles for identity genes / 11
```

This method is less dependent on absolute library scale and is used to test
whether the direction of the result is robust to scoring strategy.

## 5. Coverage and QC

- Minimum interpretable gene coverage: 70%.
- GSE176206 observed identity coverage: 11/11 genes, or 100%.
- Missing cell scores: 0.
- A score must be marked `insufficient_gene_coverage` rather than interpreted
  when coverage falls below 70%.

## 6. Biological-unit aggregation

The score is first computed per cell. For formal GSE176206 outputs, cells are
then aggregated within each unique:

```text
age_group x exact_treatment_arm x animal_label
```

The primary animal-level statistic is the median cell score. Mean and 95th
percentile values are retained as supporting summaries. This is a
`cell score -> aggregate within biological unit` pathway and is not equivalent
to scoring an animal pseudobulk expression profile.

Numeric animal labels are nested within age and exact treatment arm. Treatment
arms are therefore unpaired. Cells with `animal = unknown` are excluded from
animal-level inference and reported only in a descriptive sensitivity table.

## 7. Interpretation limits

- Higher Identity Score supports stronger expression of the frozen MSC marker
  program; it does not directly measure differentiation potency or function.
- Identity is not mathematically interchangeable with Youth or Risk.
- A treatment-associated identity decrease may reflect lineage remodeling,
  cell-state redistribution, or within-state transcriptional change.
- The small number of animals means condition comparisons are descriptive and
  should not be presented as definitive population-level inference.
