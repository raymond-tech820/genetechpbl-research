#!/usr/bin/env Rscript

local_lib <- ".Rlibs"
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(edgeR)
})

counts_path <- "data/processed/pseudobulk_filtered_counts.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata.csv"
processed_dir <- "data/processed"
qc_dir <- "outputs/qc"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

logcpm_out <- file.path(processed_dir, "pseudobulk_logcpm.rds")
logcpm_csv_out <- file.path(processed_dir, "pseudobulk_logcpm.csv")
dge_out <- file.path(processed_dir, "pseudobulk_dge_tmm.rds")
norm_factors_out <- file.path(qc_dir, "tmm_normalization_factors.csv")
verification_out <- file.path(qc_dir, "step5_normalization_verification.csv")
report_out <- file.path(qc_dir, "step5_tmm_normalization_report.md")
library_plot_out <- file.path(qc_dir, "library_sizes_before_after_normalization.png")
distribution_plot_out <- file.path(qc_dir, "logcpm_distributions_before_after_normalization.png")

message("Reading filtered pseudobulk counts")
counts <- readRDS(counts_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(colnames(counts), metadata$sample_id)) {
  stop("Filtered count columns do not match pseudobulk metadata sample_id order")
}

message("Running edgeR TMM normalization")
dge <- DGEList(counts = counts, samples = metadata)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 2)
raw_logcpm <- cpm(counts, log = TRUE, prior.count = 2)

raw_library_size <- dge$samples$lib.size
norm_factor <- dge$samples$norm.factors
effective_library_size <- raw_library_size * norm_factor

norm_table <- data.frame(
  sample_id = colnames(counts),
  age_group = metadata$age_group,
  sex = metadata$sex,
  raw_library_size = as.numeric(raw_library_size),
  tmm_norm_factor = as.numeric(norm_factor),
  effective_library_size = as.numeric(effective_library_size),
  log2_effective_over_raw = log2(effective_library_size / raw_library_size),
  stringsAsFactors = FALSE
)
write.csv(norm_table, norm_factors_out, row.names = FALSE)

message("Writing normalized expression")
saveRDS(logcpm, logcpm_out)
saveRDS(dge, dge_out)
logcpm_csv <- data.frame(gene = rownames(logcpm), logcpm, check.names = FALSE)
write.csv(logcpm_csv, logcpm_csv_out, row.names = FALSE)

message("Computing verification metrics")
sample_medians_raw <- apply(raw_logcpm, 2, median)
sample_medians_tmm <- apply(logcpm, 2, median)
sample_iqr_raw <- apply(raw_logcpm, 2, IQR)
sample_iqr_tmm <- apply(logcpm, 2, IQR)

pca_raw <- prcomp(t(raw_logcpm), center = TRUE, scale. = FALSE)
pca_tmm <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
pc1_raw_cor_library <- suppressWarnings(cor(pca_raw$x[, 1], raw_library_size, method = "spearman"))
pc1_tmm_cor_library <- suppressWarnings(cor(pca_tmm$x[, 1], raw_library_size, method = "spearman"))
pc2_raw_cor_library <- suppressWarnings(cor(pca_raw$x[, 2], raw_library_size, method = "spearman"))
pc2_tmm_cor_library <- suppressWarnings(cor(pca_tmm$x[, 2], raw_library_size, method = "spearman"))
median_raw_cor_library <- suppressWarnings(cor(sample_medians_raw, raw_library_size, method = "spearman"))
median_tmm_cor_library <- suppressWarnings(cor(sample_medians_tmm, raw_library_size, method = "spearman"))

verification <- data.frame(
  n_genes = nrow(logcpm),
  n_samples = ncol(logcpm),
  min_tmm_norm_factor = min(norm_factor),
  max_tmm_norm_factor = max(norm_factor),
  norm_factor_ratio_max_over_min = max(norm_factor) / min(norm_factor),
  any_extreme_norm_factor_lt_0_5_or_gt_2 = any(norm_factor < 0.5 | norm_factor > 2),
  raw_library_size_min = min(raw_library_size),
  raw_library_size_max = max(raw_library_size),
  effective_library_size_min = min(effective_library_size),
  effective_library_size_max = max(effective_library_size),
  sample_median_logcpm_range_before = diff(range(sample_medians_raw)),
  sample_median_logcpm_range_after = diff(range(sample_medians_tmm)),
  spearman_pc1_vs_library_before = pc1_raw_cor_library,
  spearman_pc1_vs_library_after = pc1_tmm_cor_library,
  spearman_pc2_vs_library_before = pc2_raw_cor_library,
  spearman_pc2_vs_library_after = pc2_tmm_cor_library,
  spearman_sample_median_logcpm_vs_library_before = median_raw_cor_library,
  spearman_sample_median_logcpm_vs_library_after = median_tmm_cor_library,
  sex_column_available_for_downstream_adjustment = "sex" %in% colnames(metadata),
  stringsAsFactors = FALSE
)
write.csv(verification, verification_out, row.names = FALSE)

sample_distribution <- data.frame(
  sample_id = colnames(logcpm),
  age_group = metadata$age_group,
  sex = metadata$sex,
  raw_median_logcpm = sample_medians_raw,
  tmm_median_logcpm = sample_medians_tmm,
  raw_iqr_logcpm = sample_iqr_raw,
  tmm_iqr_logcpm = sample_iqr_tmm,
  stringsAsFactors = FALSE
)
write.csv(sample_distribution, file.path(qc_dir, "step5_sample_logcpm_distribution_summary.csv"), row.names = FALSE)

message("Drawing QC plots")
png(library_plot_out, width = 1800, height = 900, res = 180)
par(mar = c(8, 5, 4, 2))
bar_positions <- barplot(
  rbind(raw_library_size, effective_library_size),
  beside = TRUE,
  col = c("#4C78A8", "#F58518"),
  names.arg = colnames(counts),
  las = 2,
  ylab = "Library size",
  main = "Library sizes before and after TMM normalization"
)
legend(
  "topright",
  legend = c("Raw filtered library size", "Effective TMM library size"),
  fill = c("#4C78A8", "#F58518"),
  bty = "n"
)
dev.off()

png(distribution_plot_out, width = 1800, height = 900, res = 180)
par(mfrow = c(1, 2), mar = c(8, 5, 4, 1))
boxplot(
  raw_logcpm,
  las = 2,
  col = "#9ECAE1",
  ylab = "log2 CPM",
  main = "Before TMM"
)
boxplot(
  logcpm,
  las = 2,
  col = "#FDBF6F",
  ylab = "log2 CPM",
  main = "After TMM"
)
dev.off()

if (any(!is.finite(logcpm))) {
  stop("Non-finite values found in normalized logCPM matrix")
}
if (!identical(colnames(logcpm), metadata$sample_id)) {
  stop("logCPM columns do not match metadata sample_id order")
}

report_lines <- c(
  "# Step 5: TMM Normalization",
  "",
  "## Inputs",
  "",
  sprintf("- Filtered pseudobulk counts: `%s`", counts_path),
  sprintf("- Pseudobulk metadata: `%s`", metadata_path),
  "",
  "## What This Step Did",
  "",
  "1. Created an edgeR `DGEList` from the filtered mouse-level pseudobulk raw counts.",
  "2. Estimated TMM normalization factors with `edgeR::calcNormFactors(method = \"TMM\")`.",
  "3. Calculated normalized log2 CPM using `edgeR::cpm(log = TRUE, prior.count = 2)`.",
  "4. Saved the normalized matrix for PCA, correlation, heatmaps, Youth Score construction, and model baselines.",
  "5. Wrote normalization-factor and expression-distribution QC summaries.",
  "",
  "## Results",
  "",
  sprintf("- Genes normalized: %s", nrow(logcpm)),
  sprintf("- Mouse-level samples: %s", ncol(logcpm)),
  sprintf("- TMM normalization factor range: %.3f to %.3f", min(norm_factor), max(norm_factor)),
  sprintf("- Max/min normalization factor ratio: %.3f", max(norm_factor) / min(norm_factor)),
  sprintf("- Any extreme normalization factor (<0.5 or >2): %s", any(norm_factor < 0.5 | norm_factor > 2)),
  sprintf("- Sample median logCPM range before TMM: %.3f", diff(range(sample_medians_raw))),
  sprintf("- Sample median logCPM range after TMM: %.3f", diff(range(sample_medians_tmm))),
  sprintf("- Spearman PC1 vs raw library size before TMM: %.3f", pc1_raw_cor_library),
  sprintf("- Spearman PC1 vs raw library size after TMM: %.3f", pc1_tmm_cor_library),
  sprintf("- Spearman sample median logCPM vs library size before TMM: %.3f", median_raw_cor_library),
  sprintf("- Spearman sample median logCPM vs library size after TMM: %.3f", median_tmm_cor_library),
  "",
  "## Sex Adjustment Note",
  "",
  "Normalization itself is unsupervised and does not include age or sex labels. The pseudobulk metadata still contains `sex`; downstream DE/modeling should include sex adjustment where the design matrix is estimable, while reporting the partial age-sex confounding noted in Step 2.",
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", logcpm_out),
  sprintf("- `%s`", logcpm_csv_out),
  sprintf("- `%s`", dge_out),
  sprintf("- `%s`", norm_factors_out),
  sprintf("- `%s`", library_plot_out),
  sprintf("- `%s`", distribution_plot_out),
  sprintf("- `%s`", verification_out)
)
writeLines(report_lines, report_out)

message("Step 5 complete")
message(sprintf("Normalized logCPM matrix: %s genes x %s samples", nrow(logcpm), ncol(logcpm)))
message(sprintf("TMM factor range: %.3f to %.3f", min(norm_factor), max(norm_factor)))
