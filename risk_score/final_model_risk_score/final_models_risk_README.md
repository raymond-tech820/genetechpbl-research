# Final Risk Score Models

This directory contains the finalized risk score definitions and the reproducible scripts used to calculate the single-cell and donor-level risk metrics for the final integration pipeline.

## Contents

### 1. `finalized_risk_gene_list.csv`
The definitive list of genes used for calculating the risk scores. The axes have been refined through rigorous biological validation:
*   **Pluripotency Risk:** Focuses exclusively on 10 strict pluripotency markers (e.g., Nanog, Pou5f1, Sox2) to accurately capture dangerous dedifferentiation events and the loss of somatic identity. Previous versions included mesenchymal identity markers (e.g., EMT), which incorrectly flagged normal somatic cells as high-risk.
*   **Genomic Stress Risk:** Derived from `HALLMARK_DNA_REPAIR` and `HALLMARK_P53_PATHWAY`. Housekeeping genes were excluded to improve specificity. (The "Inflammation" axis was removed because partial reprogramming physiologically suppresses aging-associated inflammation, meaning the score measured rejuvenation efficacy rather than intervention-induced toxicity).

### 2. `risk_score_calculate.py`
The final Python script for single-cell risk scoring using `scanpy`.
*   **Input:** Normalized and log1p-transformed single-cell expression data (e.g., `GSE176206_msc_sokm.h5ad`).
*   **Methodology:** Uses `sc.tl.score_genes` with `ctrl_size=50`, `n_bins=25`, `random_state=42`, and `use_raw=False`.
*   **Trajectory Correction:** Inverts the `velocity_pseudotime` variable (`1.0 - pseudotime`) to align with the intuitive biological direction (0.0 = baseline somatic, 1.0 = reprogrammed).
*   **Output:** An annotated `.h5ad` object with `Pluripotency_Risk_Score` and `Genomic_Stress_Risk_Score`.

### 3. `risk_score_csv_summarize.py`
The aggregation script that generates the final donor-level contract for the downstream integration module.
*   **Output Format:** An CSV table (`GSE176206_msc_sokm_risk_donor_summaries.csv`) aggregating the single-cell scores at the `mouse × treatment condition` level.
*   **Canonical Key:** Employs the required identifier structure: `GSE176206_MSC|{age_group}|{exact_treatment_arm}|animal_{animal_label}`.
*   **Metrics Calculated:** Includes central tendency (Median, Mean), upper-tail risk (95th percentile), and high-risk cell fractions evaluated against matched-control thresholds (90th, 95th, 99th percentiles), alongside gene coverage and cell count QC flags.