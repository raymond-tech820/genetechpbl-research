# FACS Youth Score v3.1: M4_raw_all_gene_denominator_pi_0_90

Role: post-ablation exploratory sensitivity comparator.

This directory is a standalone minimal package for applying one frozen model. It contains no training data and performs no fitting or recalibration.

## Files

- `signature.csv`: frozen genes, modules, weights, FACS means, and FACS standard deviations.
- `calibration.csv`: frozen young and old reference centers and coverage threshold.
- `score_facs_v3_1_frozen_model.R`: base-R pseudobulk scorer.
- `synthetic_raw_counts.csv`: synthetic genes-by-samples raw-count fixture.
- `synthetic_expected_scores.csv`: expected fixture output.
- `validate_model_package.R`: numerical self-test.
- `manifest.csv`: file sizes and MD5 hashes.

## Validate

```bash
Rscript validate_model_package.R
```

## Apply

Input must be a genes-by-samples matrix of nonnegative raw counts. Include all available genes because column sums define the frozen raw all-gene CPM denominator. Do not provide CPM, logCPM, TMM values, z-scores, or a signature-only matrix for real applications.

```r
source("score_facs_v3_1_frozen_model.R")
frame <- read.csv("path/to/raw_pseudobulk_counts.csv", check.names = FALSE)
counts <- as.matrix(frame[, -1])
storage.mode(counts) <- "numeric"
rownames(counts) <- frame[[1]]
result <- score_facs_v3_1_frozen_model(
  counts = counts,
  signature_csv = "signature.csv",
  calibration_csv = "calibration.csv",
  model = "M4_raw_all_gene_denominator_pi_0_90"
)
result$scores
result$coverage
result$missing_genes
```

Higher calibrated values indicate a more young-like FACS-relative transcriptional state. Values outside 0-1 are allowed. The score is not a probability, chronological age estimate, universal aging clock, or proof of technical independence.

Frozen signature: 29 genes (16 young-high; 13 old-high).
Minimum total and per-module weighted coverage: 0.80.
