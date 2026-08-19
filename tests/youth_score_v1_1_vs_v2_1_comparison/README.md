# FACS v1.1 vs v2.1 Model Comparison

This folder contains the frozen-model comparison between:

- Kaile Zhu: `facs_msc_youth_score_v1_1_cleaned_limb`
- Zihan Zhou: `facs_msc_youth_score_v2_1`

The comparison was run on two datasets:

1. Cleaned TMS FACS Limb Muscle MSC raw data.
2. Local TMS Droplet Limb Muscle MSC raw data.

No model was retrained, retuned, or modified. Both models were applied as frozen scorers.

## Main Report

Start here:

```text
FACS_and_Droplet_v1_1_vs_v2_1_combined_comparison_report.md
```

This combined report separates:

- internal validation evidence;
- frozen full-model same-input application;
- Droplet cross-assay sensitivity;
- paired donor bootstrap uncertainty;
- concordance decomposition;
- technical-variable audits;
- FACS gene-signature overlap.

## Subfolders

```text
facs/
droplet/
```

`facs/` contains the FACS same-input application and internal-validation checks.

`droplet/` contains the Droplet cross-assay sensitivity analysis.

## Key Interpretation Boundary

FACS same-input frozen scoring is useful for parser correctness and describing final model behavior, but it is not independent validation. Internal validation should be read from:

- Kaile v1.1 packaged donor OOF CV;
- Zihan v2.1 nested LOMO validation.

The conservative conclusion is:

Both models form directionally consistent young-old state separation on cleaned TMS FACS and local TMS Droplet donors. Kaile v1.1 shows more consistent descriptive age ordering, especially within old donors, but paired donor bootstrap does not establish statistically supported overall superiority over Zihan v2.1.

Neither model has demonstrated technical independence or external validation from these analyses.

## Important Output Tables

```text
internal_validation_summary.csv
paired_donor_bootstrap_summary.csv
model_concordance_decomposition.csv
facs/FACS_v1_1_vs_v2_1_formal_comparison_report.md
droplet/Droplet_v1_1_vs_v2_1_formal_comparison_report.md
```

`paired_donor_bootstrap_raw_iterations.csv` contains the bootstrap iterations used to compute the uncertainty intervals.
