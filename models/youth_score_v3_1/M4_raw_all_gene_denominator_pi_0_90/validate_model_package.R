#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub('^--file=', '', args[grep('^--file=', args)])
base <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
source(file.path(base, 'score_facs_v3_1_frozen_model.R'))
frame <- read.csv(file.path(base, 'synthetic_raw_counts.csv'), check.names = FALSE)
counts <- as.matrix(frame[, -1, drop = FALSE])
storage.mode(counts) <- 'numeric'
rownames(counts) <- frame[[1]]
model <- 'M4_raw_all_gene_denominator_pi_0_90'
result <- score_facs_v3_1_frozen_model(
  counts = counts,
  signature_csv = file.path(base, 'signature.csv'),
  calibration_csv = file.path(base, 'calibration.csv'),
  model = model
)$scores
expected <- read.csv(file.path(base, 'synthetic_expected_scores.csv'), check.names = FALSE)
columns <- c('raw_score', 'calibrated_score', 'young_module_score', 'old_module_score')
difference <- max(abs(as.matrix(result[, columns]) - as.matrix(expected[, columns])))
stopifnot(is.finite(difference), difference <= 1e-12, all(result$coverage_pass))
cat(sprintf('PASS: %s; max absolute difference %.3g\n', model, difference))
