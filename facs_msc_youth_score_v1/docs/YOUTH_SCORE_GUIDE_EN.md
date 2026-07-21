# TMS Youth Score: Reproducible Training and Validation Guide

## 1. Purpose

This project builds a research-grade Youth Score from mouse single-cell RNA-seq. The score is a calibrated probability between 0 and 1:

- `1.0`: expression is more similar to the TMS 3-month reference state.
- `0.0`: expression is more similar to the TMS 18–24-month reference state.

The Transformer also predicts chronological age as an auxiliary task. The age output constrains the representation but does not replace the Youth Score.

The score is not a Risk Score, causal measurement of rejuvenation, clinical age, or treatment recommendation. A high Youth Score does not establish safety or preservation of cell identity.

## 2. Dataset Design

### 2.1 Primary model

- Dataset: TMS FACS.
- Tissue: `SCAT`.
- Cell type: `mesenchymal stem cell of adipose`.
- Young: 3 months, 737 cells from 7 mice.
- Old: 18 and 24 months, 833 cells from 8 mice.

SCAT is the primary reference because its young/old cell counts and biological replicates are relatively balanced, and an adipogenic partial-reprogramming dataset is publicly available.

### 2.2 Secondary model

- Dataset: TMS FACS.
- Tissue: `Limb_Muscle`.
- Cell type: `mesenchymal stem cell`.
- Young: 3 months, 264 cells from 6 mice.
- Old: 18 and 24 months, 671 cells from 8 mice.

The Limb model provides a closer cell-type match to the external MSC experiment. It remains secondary because it has fewer cells and substantial age/sex and age/library-quality confounding.

### 2.3 Sensitivity dataset

The TMS Droplet Limb Muscle MSC set contains 13,037 cells spanning 1, 3, 18, 21, 24, and 30 months. It is not used for training because the main comparison has only 2 young versus 10 old mice, and age is strongly confounded with sex. It is used only to test cross-modality transfer and score ordering across ages.

### 2.4 External validation

The pipeline downloads two processed GSE176206 files:

- `GSE176206_adipo_sokm.h5ad.gz` for the SCAT model.
- `GSE176206_msc_sokm.h5ad.gz` for the Limb model.

External data are never used for feature selection, model fitting, early stopping, or calibration. The expected checks are:

1. Young control cells score above aged control cells.
2. Aged SOKM cells move toward the young score relative to aged controls.

Results are aggregated by biological replicate whenever replicate metadata are available. Cell-level p-values are not treated as biological replication.

Source references:

- [TMS CELLxGENE explorer](https://cellxgene.cziscience.com/e/f16a8f4d-bc97-43c5-a2f6-bbda952e4c5c.cxg/)
- [BPCells Python `DirMatrix` documentation](https://bnprks.github.io/BPCells/python/generated/bpcells.experimental.DirMatrix.html)
- [GSE176206 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176206)
- [PyTorch 2.13 release information](https://pytorch.org/blog/pytorch-2-13-release-blog/)

## 3. Why Droplet and FACS Are Separate

Droplet and FACS/Smart-seq2 differ in capture chemistry, sequencing depth, dropout structure, cell selection, gene coverage, and metadata. They are not pooled into a single training matrix. The FACS models are trained independently, while Droplet Limb is a sensitivity domain.

The tokenizer uses within-cell expression rank plus binned normalized expression. This reduces, but does not eliminate, platform differences.

## 4. Data Preparation

The local compressed BPCells archives are streamed to `data/extracted/`. The original files are not modified.

For each target population the pipeline verifies:

- Matrix dimensions.
- Exact equality and ordering of BPCells column names and metadata `index`.
- Declared tissue, cell type, age, cell, and mouse counts.
- Non-zero library size and at least 500 matrix-derived detected genes.

The prepared bundle is stored under `data/processed/youth_score/<dataset_id>/`:

- `counts.npz`: cell-by-gene sparse raw counts.
- `cells.parquet`: labels, mouse IDs, age, sex, QC, and source positions.
- `genes.parquet`: genes and primary exclusion flags.
- `donor_folds.csv`: fixed outer-fold assignments.
- `manifest.json`: source/output hashes and runtime details.

Matrix-derived counts are used because the Droplet 3-month Limb metadata have missing `n_genes` values.

## 5. Feature Handling and Confound Controls

For each outer fold, feature selection is repeated using training mice only:

1. Library-size normalization to 10,000.
2. `log1p` transformation.
3. Sparse variance calculation.
4. Selection of up to 4,096 highly variable genes.

Primary features exclude mitochondrial, ribosomal, hemoglobin, known Y-linked, `Xist`, and `Tsix` genes. These exclusions reduce obvious quality and sex shortcuts. They do not remove all sex effects.

No batch correction is applied. In TMS, age, sequencing plate, library complexity, and sex are partly confounded. Correcting a batch variable that is nearly identical to age can remove genuine aging signal or create an artificial one. Instead, the pipeline reports:

- A technical-only classifier using library size, detected genes, and sex.
- Male-only score direction.
- A 3-month versus 18-month sensitivity contrast.
- External cross-study validation.

If the selected expression model is within 0.02 ROC-AUC of the technical-only model, or a sensitivity direction reverses, the result is labeled `confound_limited`.

## 6. Transformer Architecture

Each cell is converted into a sequence containing its 256 highest-expression genes among the fold-specific vocabulary.

Each token combines:

- Learned gene identity embedding.
- 64-level expression-bin embedding.
- Within-cell expression-rank embedding.

A learned `[CLS]` token is prepended. The encoder has:

- 2 pre-normalized Transformer layers.
- Hidden width 128.
- 8 attention heads.
- Feed-forward width 512.
- GELU activation.
- Dropout 0.2.
- Random 15% gene-token dropout during training.

The two heads predict:

- Young/old logit.
- Age in months, scaled by 30 during training.

The loss is donor-balanced binary cross-entropy plus `0.25 × Huber age loss`. Every age class receives equal total weight, every mouse within a class receives equal total weight, and cells within a mouse share that mouse's weight.

## 7. Cross-validation and Calibration

The primary and secondary models use five donor-stratified outer folds. A mouse can appear in only one of train, validation, or test for a given run.

Within each outer fold:

- One young and one old training-eligible mouse are reserved for validation.
- HVGs, expression bins, model fitting, early stopping, and calibration exclude the test mice.
- Three Transformer seeds are trained and averaged.
- Temperature scaling is fitted on validation mice only.
- The binary decision threshold is selected from validation-donor mean scores by maximizing balanced accuracy; ties prefer the threshold closest to 0.5.

Out-of-fold cell probabilities are aggregated to mouse-level means. The primary metric is donor-level ROC-AUC. PR-AUC, balanced accuracy, F1, Brier score, age correlation, and a 2,000-replicate donor bootstrap interval are also reported.

## 8. Baselines and Official Model Selection

Three comparison models are produced:

- Donor-pseudobulk young/old gene signature.
- Elastic-net logistic regression.
- Technical-only diagnostic classifier.

The technical-only model cannot become the official Youth Score. Among the signature, elastic net, and Transformer:

1. Select the highest donor-level out-of-fold ROC-AUC.
2. If models are within 0.02 AUC, select the lower Brier score.
3. If still within 0.02 Brier, prefer the simpler model.

The Transformer checkpoints and predictions are retained even when a baseline wins.

## 9. Installation

From PowerShell in the project root:

```powershell
.\scripts\setup.ps1
```

For a CPU-only environment:

```powershell
.\scripts\setup.ps1 -CpuOnly
```

The setup creates `.venv`, installs pinned packages, installs the project in editable mode, and runs the environment doctor.

## 10. Commands

Check the environment:

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli doctor
```

Prepare all local datasets:

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli prepare-tms --config configs/scat_primary.yaml --config configs/limb_secondary.yaml --config configs/limb_droplet_sensitivity.yaml
```

Train and evaluate one model:

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli train --config configs/scat_primary.yaml --device cuda
.\.venv\Scripts\python.exe -m youth_score.cli evaluate --config configs/scat_primary.yaml
```

Run the all-gene feature-mask sensitivity:

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli train --config configs/scat_primary.yaml --feature-mode all --device cuda
```

Download and validate GSE176206:

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli download-external
.\.venv\Scripts\python.exe -m youth_score.cli validate-external --device cuda
```

Run everything with resumable stages:

```powershell
.\scripts\run_pipeline.ps1 -Device cuda
```

Score a prepared bundle or an h5ad file:

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli score --model-dir outputs/youth_score/scat_facs --input path/to/input.h5ad --output scores.parquet --device cuda
```

## 11. Outputs

Each trained dataset writes to `outputs/youth_score/<dataset_id>/`:

- `oof_cell_scores.parquet`.
- `oof_donor_scores.csv`.
- `model_metrics.csv`.
- `selection.json`.
- `training_report.md`.
- `training_report.json`.
- `figures/`.
- `all_genes_sensitivity/`: complete parallel OOF results with the primary exclusion mask disabled.
- `folds/fold_*/` with tokenizer state, split manifest, baseline models, calibration, Transformer checkpoints, and histories.

The public score output contains:

- `cell_id`.
- `youth_score`.
- `predicted_age_months` from the Transformer ensemble.
- `model_id` used for the Youth Score.
- `gene_overlap`.
- `qc_status`.
- Fold-specific `decision_threshold` values are retained in OOF outputs; they affect balanced accuracy/F1, not the continuous Youth Score.

## 12. Interpreting Evidence Status

- `internally_supported`: the selected model separates held-out TMS donors and does not trigger the predefined confound warning.
- `confound_limited`: technical predictors are nearly as strong or a sex/age sensitivity direction reverses.
- External support is reported separately and requires the young-control and aged-control sanity contrast before interpreting SOKM movement.

A successful software run does not guarantee biological validation. A null result, baseline victory, or confound warning is a valid and important outcome.

## 13. Completed Reference Run (2026-07-17)

The checked-in workflow was executed on the complete target cohorts, not only on synthetic data.

### 13.1 Internal donor-aware validation

| Model domain | Official model | Donor ROC-AUC | Brier | Transformer ROC-AUC | Technical-only ROC-AUC | Status |
|---|---:|---:|---:|---:|---:|---|
| FACS SCAT adipose MSC | Elastic net | 0.9464 | 0.0883 | 0.8571 | 0.6964 | `internally_supported` |
| FACS Limb Muscle MSC | Gene signature | 1.0000 | 0.0467 | 0.8750 | 0.6875 | `internally_supported` |

The Transformer was therefore retained as a fully trained comparison and auxiliary age predictor, but it was not forced to become the official Youth Score. The Droplet Limb sensitivity score had donor-level Spearman correlation -0.8564 with chronological age, so its overall age direction passed; this remains a confounded cross-modality sensitivity result.

The all-gene sensitivity runs retained the excluded mitochondrial, ribosomal, hemoglobin, Y-linked, `Xist`, and `Tsix` genes. They selected the same official model in both domains: SCAT elastic net remained at AUC 0.9464, and Limb gene signature remained at AUC 1.0000. The SCAT Transformer decreased from AUC 0.8571 to 0.6964, while the Limb Transformer remained at 0.8750. Thus the official model choices were robust to this feature-mask sensitivity, although the SCAT Transformer representation was not.

### 13.2 GSE176206 preflight and validation

| External domain | Minimum gene overlap | TMS reference-expression Spearman | Young control - aged control | Aged SOKM - aged control | Replication |
|---|---:|---:|---:|---:|---|
| Adipo / SCAT model | 0.9431 | 0.5710 | +0.0498 | -0.0109 | One pooled library per condition; no donor CI |
| MSC / Limb model | 0.9421 | 0.5319 | +0.0298 (95% CI +0.0135 to +0.0453) | -0.1185 (95% CI -0.1612 to -0.0551) | Three identifiable animals per condition |

Both external control contrasts passed the basic age-direction check. Neither SOKM contrast moved aged cells toward a higher Youth Score; the MSC result moved significantly in the opposite direction under this model. This run therefore supports cross-study recognition of the young-versus-aged control direction, but it does **not** validate SOKM-associated rejuvenation by this Youth Score.

The Adipo h5ad does not expose individual-animal identifiers, so its condition-level result is descriptive and labeled `replication_limited`. The MSC h5ad uses `animal` as the biological unit; donor-unknown cells are scored at cell level but excluded from donor summaries. Only `Tg+/Dox+` is mapped to SOKM; `Tg+/Dox-`, `Tg-/Dox+`, and `NegCtrl` are mapped to control.

## 14. Troubleshooting

- If BPCells import fails, rerun `scripts/setup.ps1`; the project requires the pinned Windows CPython 3.12 wheel.
- If CUDA is unavailable, run `doctor`, confirm the NVIDIA driver, or use `--device cpu` for testing.
- If a download is interrupted, rerun `download-external`; `.part` files resume using HTTP byte ranges.
- If external gene overlap is below 70%, inspect `metadata_schema.json` and the h5ad gene-symbol column. Scoring stops rather than silently producing an invalid result.
- If `confound_limited` is reported, do not remove the warning. Review the technical-only model and external control contrast before using the score.
