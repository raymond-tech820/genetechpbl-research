# Risk Score Reports & Validations

This directory contains the final outputs, summary tables, and validation plots supporting the construct validity of the finalized Risk Score modules.

## Contents

### 1. `GSE176206_msc_sokm_risk_donor_summaries.csv`
The definitive donor-level aggregation table ready for cross-module integration.
*   **Canonical Key:** `GSE176206_MSC|{age_group}|{exact_treatment_arm}|animal_{animal_label}`
*   Contains QC flags (cell counts, gene coverage) and primary risk metrics (Median, Mean, p95, and high-risk fractions at p90, p95, and p99) for both the Pluripotency and Genomic Stress axes across all biological replicates.

### 2. `Finalized_risk_score_roux.png`
Scatter plots with LOWESS regression lines depicting the finalized Risk Scores (Pluripotency and Genomic Stress) across the corrected reprogramming trajectory.
*   **Data Source:** Roux et al., 2022 (*Cell Systems* 13, 574–587) - "Diverse partial reprogramming strategies restore youthful gene expression and transiently suppress
cell identity" (GEO: GSE176206).
*   **Key Insight:** Demonstrates the correct temporal dynamics where the Pluripotency risk spikes at the distal end of the reprogramming trajectory. (The Inflammation axis was intentionally omitted here, as it measures rejuvenation efficacy rather than toxicity).

### 3. `pluripotency_Positive_Control_Validation.png`
Positive control validation for the Pluripotency Risk Score.
*   **Data Source:** Velychko et al., 2019 (*Cell Stem Cell* 25, 737–753) - "Excluding Oct4 from Yamanaka Cocktail Unleashes the Developmental Potential of iPSCs" (GEO: GSE137001).
*   **Key Insight:** Confirms construct validity by showing a massive, accurate spike in the score exclusively in fully dedifferentiated states (iPSCs/ESCs) compared to somatic MEFs and reprogramming intermediates.

### 4. `Genomic_Stress_positive_control_PDGFRA.png`
Positive control validation for the Genomic Stress Risk Score.
*   **Data Source:** Naipauer et al., 2019 (*PLOS Pathogens*) - "PDGFRA defines the mesenchymal stem cell Kaposi’s sarcoma progenitors by enabling KSHV oncogenesis in an angiogenic environment" (GEO: GSE141866).
*   **Key Insight:** Validates the axis by demonstrating a significant score elevation in actual *in vivo* sarcoma tumor samples compared to normal MSCs and pre-cancerous intermediates.