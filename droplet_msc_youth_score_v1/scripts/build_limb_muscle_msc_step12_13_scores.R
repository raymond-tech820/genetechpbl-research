#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stats)
})

young_rank_path <- "outputs/signature/step10_young_high_ranked_candidates.csv"
old_rank_path <- "outputs/signature/step10_old_high_ranked_candidates.csv"
logcpm_path <- "data/processed/pseudobulk_logcpm.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv"
scores_dir <- "outputs/scores"
models_dir <- "models"

dir.create(scores_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

w_max <- 3
versions <- data.frame(
  version = c("Small", "Medium", "Large"),
  genes_per_direction = c(20L, 50L, 100L),
  stringsAsFactors = FALSE
)

required_rank_cols <- c("gene", "direction", "adjusted_logFC", "age_rho", "r_g", "q_g")

read_rank <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required_rank_cols, colnames(x))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns in %s: %s", path, paste(missing, collapse = ", ")))
  }
  x
}

weighted_mean_safe <- function(z, w) {
  ok <- is.finite(z) & is.finite(w) & abs(w) > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(w[ok] * z[ok]) / sum(abs(w[ok]))
}

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

format_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

young_rank <- read_rank(young_rank_path)
old_rank <- read_rank(old_rank_path)
logcpm <- readRDS(logcpm_path)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(metadata$sample_id %in% colnames(logcpm))) {
  missing_samples <- setdiff(metadata$sample_id, colnames(logcpm))
  stop(sprintf("Metadata samples missing from logCPM matrix: %s", paste(missing_samples, collapse = ", ")))
}

metadata <- metadata[match(colnames(logcpm), metadata$sample_id), ]
if (!all(metadata$sample_id == colnames(logcpm))) {
  stop("Metadata could not be aligned to logCPM columns.")
}

selected_signatures <- list()
score_rows <- list()
coverage_rows <- list()
calibration_rows <- list()
signature_export_rows <- list()

for (i in seq_len(nrow(versions))) {
  version_name <- versions$version[i]
  n_dir <- versions$genes_per_direction[i]

  young_sel <- head(young_rank, n_dir)
  old_sel <- head(old_rank, n_dir)
  sig <- rbind(young_sel, old_sel)
  sig$module <- ifelse(sig$direction == "young_high", "young_module", "old_module")
  sig$version <- version_name
  sig$genes_per_direction <- n_dir
  sig$weight <- pmin(abs(sig$adjusted_logFC), w_max) * sig$r_g
  sig$training_mean <- NA_real_
  sig$training_sd <- NA_real_
  sig$available_in_logcpm <- sig$gene %in% rownames(logcpm)

  for (j in seq_len(nrow(sig))) {
    gene <- sig$gene[j]
    if (gene %in% rownames(logcpm)) {
      values <- as.numeric(logcpm[gene, ])
      sig$training_mean[j] <- mean(values, na.rm = TRUE)
      sig$training_sd[j] <- sd(values, na.rm = TRUE)
    }
  }

  sig$usable_for_score <- sig$available_in_logcpm & is.finite(sig$training_sd) & sig$training_sd > 0 &
    is.finite(sig$weight) & sig$weight > 0
  selected_signatures[[version_name]] <- sig
  signature_export_rows[[version_name]] <- sig

  young_sig <- sig[sig$module == "young_module" & sig$usable_for_score, ]
  old_sig <- sig[sig$module == "old_module" & sig$usable_for_score, ]
  if (nrow(young_sig) == 0 || nrow(old_sig) == 0) {
    stop(sprintf("No usable genes for %s in one or both modules.", version_name))
  }

  sample_scores <- metadata
  sample_scores$version <- version_name
  sample_scores$score_young_module_raw <- NA_real_
  sample_scores$score_old_module_raw <- NA_real_
  sample_scores$score_raw <- NA_real_

  for (sample_id in metadata$sample_id) {
    young_z <- (as.numeric(logcpm[young_sig$gene, sample_id]) - young_sig$training_mean) / young_sig$training_sd
    old_z <- (as.numeric(logcpm[old_sig$gene, sample_id]) - old_sig$training_mean) / old_sig$training_sd
    row_index <- which(sample_scores$sample_id == sample_id)
    sample_scores$score_young_module_raw[row_index] <- weighted_mean_safe(young_z, young_sig$weight)
    sample_scores$score_old_module_raw[row_index] <- weighted_mean_safe(old_z, old_sig$weight)
    sample_scores$score_raw[row_index] <- sample_scores$score_young_module_raw[row_index] -
      sample_scores$score_old_module_raw[row_index]
  }

  young_center <- median(sample_scores$score_raw[sample_scores$age_group == "Young"], na.rm = TRUE)
  old_center <- median(sample_scores$score_raw[sample_scores$age_group == "Old"], na.rm = TRUE)
  denom <- young_center - old_center
  if (!is.finite(denom) || abs(denom) < .Machine$double.eps) {
    stop(sprintf("Cannot calibrate %s: young and old centers are identical or invalid.", version_name))
  }

  sample_scores$young_reference_center <- young_center
  sample_scores$old_reference_center <- old_center
  sample_scores$youth_score_raw_calibrated <- (sample_scores$score_raw - old_center) / denom
  sample_scores$youth_score_clipped_0_1 <- clip01(sample_scores$youth_score_raw_calibrated)

  score_rows[[version_name]] <- sample_scores

  total_weight <- sum(abs(sig$weight[sig$usable_for_score]))
  young_weight <- sum(abs(sig$weight[sig$usable_for_score & sig$module == "young_module"]))
  old_weight <- sum(abs(sig$weight[sig$usable_for_score & sig$module == "old_module"]))

  for (sample_id in metadata$sample_id) {
    available <- sig$usable_for_score & !is.na(logcpm[sig$gene, sample_id])
    available_young <- available & sig$module == "young_module"
    available_old <- available & sig$module == "old_module"
    coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
      version = version_name,
      sample_id = sample_id,
      signature_genes_total = nrow(sig),
      signature_genes_available = sum(available),
      total_coverage = sum(available) / nrow(sig),
      young_module_genes_total = sum(sig$module == "young_module"),
      young_module_genes_available = sum(available_young),
      young_module_coverage = sum(available_young) / sum(sig$module == "young_module"),
      old_module_genes_total = sum(sig$module == "old_module"),
      old_module_genes_available = sum(available_old),
      old_module_coverage = sum(available_old) / sum(sig$module == "old_module"),
      weighted_coverage = sum(abs(sig$weight[available])) / total_weight,
      young_weighted_coverage = sum(abs(sig$weight[available_young])) / young_weight,
      old_weighted_coverage = sum(abs(sig$weight[available_old])) / old_weight,
      stringsAsFactors = FALSE
    )
  }

  calibration_rows[[version_name]] <- data.frame(
    version = version_name,
    genes_per_direction = n_dir,
    young_high_genes = nrow(young_sig),
    old_high_genes = nrow(old_sig),
    young_reference_center = young_center,
    old_reference_center = old_center,
    calibration_denominator = denom,
    score_weight_formula = "min(abs(adjusted_logFC), 3) * r_g",
    ranking_formula_from_step10 = "|adjusted_logFC| * |age_rho| * r_g",
    score_interpretation = "full-data apparent training score; not leakage-safe validation",
    stringsAsFactors = FALSE
  )
}

signature_versions <- do.call(rbind, signature_export_rows)
candidate_scores <- do.call(rbind, score_rows)
gene_coverage <- do.call(rbind, coverage_rows)
calibration <- do.call(rbind, calibration_rows)

signature_cols <- c(
  "version", "genes_per_direction", "gene", "module", "direction",
  "adjusted_logFC", "age_rho", "r_g", "q_g", "weight",
  "training_mean", "training_sd", "available_in_logcpm", "usable_for_score"
)
extra_cols <- setdiff(colnames(signature_versions), signature_cols)
signature_versions <- signature_versions[, c(signature_cols, extra_cols), drop = FALSE]

score_cols <- c(
  "version", "sample_id", "mouse_id", "age_months", "age_group", "sex",
  "cell_count", "raw_library_size", "effective_library_size",
  "score_young_module_raw", "score_old_module_raw", "score_raw",
  "young_reference_center", "old_reference_center",
  "youth_score_raw_calibrated", "youth_score_clipped_0_1"
)
candidate_scores <- candidate_scores[, score_cols, drop = FALSE]

write.csv(signature_versions, file.path(scores_dir, "step12_candidate_signature_versions.csv"), row.names = FALSE)
write.csv(candidate_scores, file.path(scores_dir, "step12_13_candidate_scores.csv"), row.names = FALSE)
write.csv(gene_coverage, file.path(scores_dir, "step13_gene_coverage.csv"), row.names = FALSE)
write.csv(calibration, file.path(scores_dir, "step13_training_calibration.csv"), row.names = FALSE)

medium_signature <- signature_versions[signature_versions$version == "Medium", ]
write.csv(
  medium_signature,
  file.path(models_dir, "limb_msc_general_youth_score_v1_medium_signature_candidates_step12.csv"),
  row.names = FALSE
)

version_summary <- aggregate(
  youth_score_raw_calibrated ~ version + age_group,
  candidate_scores,
  function(x) median(x, na.rm = TRUE)
)
names(version_summary)[names(version_summary) == "youth_score_raw_calibrated"] <- "median_youth_score_raw_calibrated"

score_range <- aggregate(
  youth_score_raw_calibrated ~ version,
  candidate_scores,
  function(x) paste(format_num(min(x, na.rm = TRUE)), format_num(max(x, na.rm = TRUE)), sep = " to ")
)
names(score_range)[names(score_range) == "youth_score_raw_calibrated"] <- "raw_calibrated_range"

coverage_summary <- aggregate(
  cbind(total_coverage, young_module_coverage, old_module_coverage, weighted_coverage) ~ version,
  gene_coverage,
  min
)

report_lines <- c(
  "# Step 12-13: Candidate Youth Score Construction and Training Calibration",
  "",
  "## Inputs",
  "",
  sprintf("- Young-high ranking: `%s`", young_rank_path),
  sprintf("- Old-high ranking: `%s`", old_rank_path),
  sprintf("- TMM logCPM matrix: `%s`", logcpm_path),
  sprintf("- Mouse metadata: `%s`", metadata_path),
  "",
  "## What This Step Did",
  "",
  "1. Built three candidate signatures from the Step 10 compressed rankings: Small = 20+20, Medium = 50+50, Large = 100+100.",
  "2. Computed final scoring weights as `min(abs(adjusted_logFC), 3) * r_g`.",
  "3. Did not multiply `abs(age_rho)` into the final scoring weight, because Step 10 already used it in `q_g` for ranking.",
  "4. Standardized each signature gene across the 12 mouse pseudobulk samples with full-training `mu_g` and `sd_g`.",
  "5. Computed raw module scores as weighted averages of z-scored expression, then `S = S_young - S_old`.",
  "6. Calibrated raw scores using training medians: `Y_raw = (S - M_O) / (M_Y - M_O)` and also exported `Y_clipped_0_1`.",
  "7. Exported total, young-module, old-module, and weighted gene coverage.",
  "",
  "## Important Interpretation",
  "",
  "These are full-data apparent training scores. They are useful for checking score construction, calibration, and coverage, but they are not leakage-safe validation. Step 14 must redo feature selection inside each held-out-mouse fold before using LOMO performance as evidence.",
  "",
  "## Calibration Summary",
  "",
  paste(capture.output(print(calibration, row.names = FALSE)), collapse = "\n"),
  "",
  "## Median Calibrated Score by Age Group",
  "",
  paste(capture.output(print(version_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Calibrated Score Range",
  "",
  paste(capture.output(print(score_range, row.names = FALSE)), collapse = "\n"),
  "",
  "## Minimum Gene Coverage Across Samples",
  "",
  paste(capture.output(print(coverage_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Outputs",
  "",
  "- `outputs/scores/step12_candidate_signature_versions.csv`",
  "- `outputs/scores/step12_13_candidate_scores.csv`",
  "- `outputs/scores/step13_gene_coverage.csv`",
  "- `outputs/scores/step13_training_calibration.csv`",
  "- `models/limb_msc_general_youth_score_v1_medium_signature_candidates_step12.csv`",
  "",
  "## Next Step",
  "",
  "Proceed to Step 14 with nested leave-one-mouse-out validation. Do not reuse these full-data signatures to score held-out mice in that validation."
)

writeLines(report_lines, file.path(scores_dir, "step12_13_candidate_score_report.md"))

cat("Step 12-13 complete\n")
print(calibration, row.names = FALSE)
cat("\nMedian calibrated score by age group:\n")
print(version_summary, row.names = FALSE)
cat("\nMinimum gene coverage across samples:\n")
print(coverage_summary, row.names = FALSE)
