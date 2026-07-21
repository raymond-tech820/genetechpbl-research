# FACS MSC Youth Score v1

This deliverable contains the reproducible FACS-trained Youth Score models used in the Genetech PBL project. It provides two mouse single-cell RNA-seq reference models from Tabula Muris Senis (TMS):

| Reference cohort | Selected Youth Score model | Cells | Donors | Donor-level OOF ROC-AUC |
| --- | --- | ---: | ---: | ---: |
| SCAT adipose MSC | Elastic Net | 1,570 | 15 | 0.946 |
| Limb Muscle MSC | Gene-signature score | 935 | 14 | 1.000 |

Both models distinguish 3-month mice from 18- and 24-month mice. Scores are calibrated probabilities of similarity to the young TMS reference, not clinical ages, safety scores, or causal measures of rejuvenation.

## Quick start

The scorer expects a cell-by-gene `.h5ad` object containing **raw, non-negative count data** in `adata.raw.X` (preferred) or `adata.X`. Do not provide log-normalized counts, CPM/TPM, or scaled data. The scorer performs per-cell library-size normalization to 10,000 and `log1p` internally.

```powershell
cd research_pipeline
..\\scripts\\setup.ps1

# Limb Muscle MSC selected model
youth-score score --model-dir ..\\models\\limb_facs --input your_counts.h5ad --output limb_scores.csv

# SCAT adipose MSC selected model
youth-score score --model-dir ..\\models\\scat_facs --input your_counts.h5ad --output scat_scores.csv
```

The output contains `cell_id`, `youth_score`, `predicted_age_months`, `model_id`, `gene_overlap`, and `qc_status`. A `gene_overlap` below 0.70 fails the model's preflight requirement.

## What is included

- `models/`: all fold-specific tokenizer, calibration, baseline, and Transformer checkpoint assets needed by the original scorer.
- `research_pipeline/`: English Python implementation, configurations, and tests. Raw TMS and GSE176206 data are intentionally excluded.
- `results/`: donor-level cross-validation scores, Droplet cross-platform sensitivity results, and external GSE176206 summary tables.
- `figures/`: summary figures used for the project report.
- `docs/`: model I/O documentation and bilingual project guides.

## Key validation results

- Cross-platform sensitivity: the FACS Limb model scored TMS Droplet Limb MSC donors with Spearman rho = -0.856 across 16 mice (age 1 to 30 months).
- External GSE176206 controls: young controls scored above aged controls for both the matched Adipo and MSC datasets.
- Aged SOKM samples did not move in the young direction under this TMS-defined score. This is a negative result for this score, not evidence about every possible rejuvenation phenotype.

## Important limitations

Validation is donor-aware, but the TMS cohorts are small and age, sex, library complexity, and assay can be partially confounded. Aggregate scores by true donor or biological replicate; never treat cells as independent biological replicates. Compare scores primarily within a dataset, rather than interpreting an absolute score threshold across platforms.

See [MODEL_CARD.md](MODEL_CARD.md) and [docs/YOUTH_SCORE_MODEL_INPUT_OUTPUT.md](docs/YOUTH_SCORE_MODEL_INPUT_OUTPUT.md) before applying the models.
