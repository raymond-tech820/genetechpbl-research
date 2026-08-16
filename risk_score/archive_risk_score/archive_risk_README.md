# Archive: Evolution of the Risk Score Module

This directory serves as the reproducible record of the Risk Score module's development. It contains previous iterations of the gene lists and exploratory plots that highlight the biological reasoning behind the final model choices.

## Historical Gene Lists

1.  **`risk_genes.csv` (Version 1)**
    *   **Description:** The initial compilation mapping raw MSigDB HALLMARK pathways directly to three axes (Dedifferentiation, Inflammation, Tumor).
    *   **Issues Identified:** 
        *   **Dedifferentiation:** Included mesenchymal markers (e.g., EMT, Wnt, Hedgehog, Notch) alongside 10 pluripotency markers. This incorrectly penalized normal somatic cells.
        *   **Tumor:** Contained over 800 genes related to normal cellular proliferation (E2F targets, Myc targets, G2M checkpoint, Mitotic spindle, p53 pathway).

2.  **`risk_genes_adjusted.csv` (Version 2)**
    *   **Description:** The "Tumor" axis was renamed to "Genomic Stress" and narrowed down to `HALLMARK_DNA_REPAIR` and `HALLMARK_P53_PATHWAY` (338 genes).
    *   **Issues Resolved:** Prevented the score from falsely flagging the normal temporary cell cycle arrest caused by SOKM stress as a tumorigenic event.

3.  **`risk_genes_housekeeping_revised.csv` (Version 3)**
    *   **Description:** Filtered out genes present in the `HOUNKPE_HOUSEKEEPING_GENES` dataset from the Version 2 list to increase the specificity of the risk signals.

## Exploratory Plots and Justifications

### 1. `risk_score_dediff_compare.png`
*   **Description:** Compares the old combined "Dediff" score against scoring only the "Mesenchymal" markers and only the "Pluripotency" markers.
*   **Key Insight:** Proves that the original Dediff score was heavily skewed by mesenchymal markers, showing a false "decrease in risk" as reprogramming progressed. This justified narrowing the axis exclusively to the 10 Pluripotency markers.

### 2. `Finalized_risk_score_roux_with_inflam.png`
*   **Description:** Plots the finalized genes but includes the Inflammation axis across the corrected pseudotime trajectory.
*   **Key Insight:** Demonstrates that the Inflammation score *decreases* during reprogramming. This aligns with Roux et al.'s findings that partial reprogramming suppresses age-associated inflammation. It justified removing the Inflammation axis from the "Risk" module, as it was capturing rejuvenation efficacy rather than intervention-induced toxicity.

### 3. `risk_score_tumor_compare.png`
*   **Description:** Compares the old "Tumor" score with the new "Genomic Stress" score.
*   **Key Insight:** The old Tumor score spiked late in the trajectory, which likely reflected the physiological resumption of normal cell division rather than actual tumorigenesis. The shift to Genomic Stress resolved this artifact.

### 4. `risk_score_genomic_stress_housekeeping_compare.png`
*   **Description:** Compares the Genomic Stress score with and without housekeeping gene filtering.
*   **Key Insight:** Shows minimal difference in the final output distributions. (Pluripotency was not plotted here as it naturally contained no housekeeping genes).
