# Youth Score Model: Input and Output Specification

## 1. Model variants

This project provides two tissue-specific Youth Score ensembles:

| Model directory | Intended cell population | Official Youth Score model |
|---|---|---|
| `outputs/youth_score/scat_facs` | SCAT adipose MSC-like cells | Elastic Net |
| `outputs/youth_score/limb_facs` | Limb Muscle MSCs | Gene signature |

The models were trained on TMS FACS single-cell raw counts. They should not be treated as one universal pan-cell-type aging model.

## 2. Required input

### 2.1 Supported input formats

The command-line scorer accepts either:

1. An `.h5ad` file; or
2. A prepared Youth Score bundle directory containing `counts.npz`, `cells.parquet`, and `genes.parquet`.

For an `.h5ad` file, the scorer reads:

- Expression matrix: `adata.raw.X` when `adata.raw` exists; otherwise `adata.X`.
- Gene identifiers: the corresponding `var_names`.
- Cell identifiers: `adata.obs_names`.

The required matrix orientation is:

```text
rows    = cells
columns = genes
values  = raw, non-negative expression counts
```

Sparse matrices are preferred. Each row should represent one target cell, not a donor-level pseudobulk sample.

### 2.2 Raw counts are required

The input should contain raw UMI/read counts. Do **not** provide:

- library-size-normalized values;
- CPM, logCPM, TPM, FPKM, or RPKM;
- `log1p`-transformed expression;
- z-scores or scaled expression;
- batch-corrected or integrated expression values.

If `adata.raw` exists, verify that `adata.raw.X` truly contains raw counts. Some workflows use `.raw` to store log-normalized expression; such an object would be unsuitable and would be normalized twice.

### 2.3 Gene identifiers

Gene matching is performed by exact gene-name strings against each training-fold vocabulary.

- Mouse gene symbols matching the TMS vocabulary are recommended.
- Ensembl identifiers are not converted automatically.
- Gene-name case is not corrected automatically.
- Extra genes are ignored.
- Missing model genes are filled with zero.
- Duplicate input gene names should be collapsed before scoring.

A minimum vocabulary overlap of 70% is required for a passing QC status. Higher overlap is preferable.

### 2.4 Cell selection and QC

Input cells should already represent the intended target population:

- SCAT adipose MSC-like cells for the SCAT model; or
- Limb Muscle MSCs for the Limb model.

Recommended cell-level requirements are:

- non-zero total counts;
- at least approximately 500 detected genes, matching the training filter;
- removal of low-quality cells, doublets, and clear cell-type contaminants before scoring.

The scorer does not perform cell-type annotation, doublet detection, or batch correction.

## 3. Internal preprocessing

The user should not normalize the input. The scorer performs the required preprocessing internally and independently for every cell:

1. Align input genes to each fold-specific training vocabulary.
2. Normalize each cell to a total library size of 10,000.
3. Apply `log1p` transformation.
4. Use the fixed training-fold features and calibration parameters.
5. Average predictions across the five cross-validation folds.

For the Transformer auxiliary model, the scorer additionally:

- uses up to 4,096 fold-specific highly variable genes learned during training;
- retains the 256 highest-expression model genes per cell;
- converts normalized expression into 64 training-derived expression bins;
- encodes gene identity, expression rank, and expression bin.

No new highly variable genes, expression bins, model weights, or calibration parameters are fitted on the scoring dataset. No batch correction is applied.

## 4. Output

The output is a CSV or Parquet table, depending on the requested file extension, with one row per input cell.

| Column | Meaning |
|---|---|
| `cell_id` | Input cell identifier from `obs_names` or the prepared bundle |
| `youth_score` | Continuous score between 0 and 1, averaged across five trained folds; higher values indicate greater similarity to the TMS young-expression state |
| `predicted_age_months` | Auxiliary month-age prediction from the Transformer ensemble |
| `model_id` | Official model used for Youth Score, normally `elastic_net` for SCAT or `gene_signature` for Limb |
| `gene_overlap` | Minimum fraction of required model genes found across the five folds |
| `qc_status` | `pass` when gene overlap is at least 0.70; otherwise `insufficient_gene_overlap` |

The official `youth_score` may come from Elastic Net or Gene signature even though `predicted_age_months` always comes from the retained Transformer ensemble.

## 5. Interpretation

- A higher Youth Score means that a cell is more similar to the young TMS reference along the learned transcriptomic axis.
- The score is a research-grade relative transcriptomic score, not a clinical age, safety score, causal rejuvenation measurement, or treatment-efficacy endpoint.
- Absolute calibration may shift across assays, laboratories, and cell states. For external datasets, prioritize within-dataset rankings and contrasts instead of applying a universal 0.5 cutoff.
- Cells are not independent biological replicates. For group comparisons, first aggregate cell-level scores by true donor or experimental replicate, then compare donor-level values.
- `predicted_age_months` is an auxiliary output and should be interpreted more cautiously than the selected Youth Score.
- A passing gene-overlap QC status confirms vocabulary coverage only; it does not confirm cell identity or overall data quality.

## 6. Command-line example

From the project root:

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli score `
  --model-dir outputs/youth_score/limb_facs `
  --input path/to/limb_msc_raw_counts.h5ad `
  --output path/to/limb_msc_youth_scores.parquet `
  --device cuda
```

Use `--device cpu` when CUDA is unavailable. Use the SCAT model directory for appropriately annotated SCAT adipose MSC-like cells.

## 7. Important distinction from a pseudobulk model

This scorer was trained and implemented for single-cell, cell-by-gene raw-count input. It is not the same interface as a gene-by-sample pseudobulk logCPM model. A pseudobulk matrix should not be passed directly to this scorer without a separately validated adaptation.

