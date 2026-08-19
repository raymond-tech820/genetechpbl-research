#!/usr/bin/env Rscript

local_lib <- ".Rlibs"
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

logcpm_path <- "data/processed/pseudobulk_logcpm.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata.csv"
norm_factors_path <- "outputs/qc/tmm_normalization_factors.csv"
eda_dir <- "outputs/eda"
qc_dir <- "outputs/qc"

dir.create(eda_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

top_variable_gene_count <- 2000L

message("Reading normalized pseudobulk data")
logcpm <- readRDS(logcpm_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
norm_factors <- read.csv(norm_factors_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(colnames(logcpm), metadata$sample_id)) {
  stop("logCPM columns do not match metadata sample_id order")
}
if (!identical(metadata$sample_id, norm_factors$sample_id)) {
  stop("metadata sample_id order does not match normalization-factor table")
}

sample_info <- merge(
  metadata,
  norm_factors[, c("sample_id", "raw_library_size", "tmm_norm_factor", "effective_library_size", "log2_effective_over_raw")],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
sample_info <- sample_info[match(metadata$sample_id, sample_info$sample_id),]
sample_info$age_binary_young <- ifelse(sample_info$age_group == "Young", 1, 0)
sample_info$sex_binary_male <- ifelse(sample_info$sex == "male", 1, 0)

message("Selecting high-variance genes for PCA and heatmaps")
gene_variance <- apply(logcpm, 1, var)
top_n <- min(top_variable_gene_count, length(gene_variance))
top_variable_genes <- names(sort(gene_variance, decreasing = TRUE))[seq_len(top_n)]
write.csv(
  data.frame(gene = top_variable_genes, variance = gene_variance[top_variable_genes], row.names = NULL),
  file.path(eda_dir, "top_variable_genes_for_pca.csv"),
  row.names = FALSE
)

eda_matrix <- logcpm[top_variable_genes, , drop = FALSE]
pca <- prcomp(t(eda_matrix), center = TRUE, scale. = FALSE)
pct_var <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
pca_table <- cbind(
  sample_info,
  data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    stringsAsFactors = FALSE
  )
)
write.csv(pca_table, file.path(eda_dir, "pseudobulk_pca_scores.csv"), row.names = FALSE)

message("Computing PC1 associations")
safe_spearman <- function(x, y) {
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(estimate = unname(test$estimate), p_value = test$p.value)
}
safe_wilcox <- function(pc, group) {
  test <- suppressWarnings(wilcox.test(pc ~ group, exact = FALSE))
  c(statistic = unname(test$statistic), p_value = test$p.value)
}

pc1_age_months <- safe_spearman(pca_table$PC1, pca_table$age_months)
pc1_age_binary <- safe_wilcox(pca_table$PC1, pca_table$age_group)
pc1_sex <- safe_wilcox(pca_table$PC1, pca_table$sex)
pc1_cell_count <- safe_spearman(pca_table$PC1, pca_table$cell_count)
pc1_raw_library <- safe_spearman(pca_table$PC1, pca_table$raw_library_size)
pc1_effective_library <- safe_spearman(pca_table$PC1, pca_table$effective_library_size)

pc1_associations <- data.frame(
  variable = c(
    "age_months",
    "age_group_binary",
    "sex",
    "cell_count",
    "raw_library_size",
    "effective_library_size"
  ),
  method = c("Spearman", "Wilcoxon rank-sum", "Wilcoxon rank-sum", "Spearman", "Spearman", "Spearman"),
  statistic = c(
    pc1_age_months["estimate"],
    pc1_age_binary["statistic"],
    pc1_sex["statistic"],
    pc1_cell_count["estimate"],
    pc1_raw_library["estimate"],
    pc1_effective_library["estimate"]
  ),
  p_value = c(
    pc1_age_months["p_value"],
    pc1_age_binary["p_value"],
    pc1_sex["p_value"],
    pc1_cell_count["p_value"],
    pc1_raw_library["p_value"],
    pc1_effective_library["p_value"]
  ),
  stringsAsFactors = FALSE
)
write.csv(pc1_associations, file.path(eda_dir, "pc1_associations.csv"), row.names = FALSE)

message("Auditing extreme PC1 and library-size samples")
pca_table$abs_PC1_z <- abs(as.numeric(scale(pca_table$PC1)))
pca_table$abs_effective_library_z <- abs(as.numeric(scale(pca_table$effective_library_size)))
pca_table$abs_raw_library_z <- abs(as.numeric(scale(pca_table$raw_library_size)))
pca_table$extreme_PC1_top2 <- pca_table$sample_id %in% pca_table$sample_id[order(-pca_table$abs_PC1_z)][1:2]
pca_table$extreme_effective_library_top2 <- pca_table$sample_id %in% pca_table$sample_id[order(-pca_table$abs_effective_library_z)][1:2]
pca_table$extreme_raw_library_top2 <- pca_table$sample_id %in% pca_table$sample_id[order(-pca_table$abs_raw_library_z)][1:2]
pca_table$extreme_PC1_and_effective_library <- pca_table$extreme_PC1_top2 & pca_table$extreme_effective_library_top2
pca_table$extreme_PC1_and_raw_library <- pca_table$extreme_PC1_top2 & pca_table$extreme_raw_library_top2
extreme_audit <- pca_table[
  order(-pca_table$abs_PC1_z),
  c(
    "sample_id",
    "age_group",
    "sex",
    "cell_count",
    "raw_library_size",
    "effective_library_size",
    "tmm_norm_factor",
    "PC1",
    "PC2",
    "abs_PC1_z",
    "abs_effective_library_z",
    "abs_raw_library_z",
    "extreme_PC1_top2",
    "extreme_effective_library_top2",
    "extreme_raw_library_top2",
    "extreme_PC1_and_effective_library",
    "extreme_PC1_and_raw_library"
  )
]
write.csv(extreme_audit, file.path(eda_dir, "extreme_pc1_library_audit.csv"), row.names = FALSE)

library_audit <- pca_table[
  ,
  c(
    "sample_id",
    "age_group",
    "sex",
    "raw_library_size",
    "effective_library_size",
    "tmm_norm_factor",
    "log2_effective_over_raw",
    "PC1",
    "PC2"
  )
]
library_audit$effective_minus_raw <- library_audit$effective_library_size - library_audit$raw_library_size
library_audit$effective_over_raw <- library_audit$effective_library_size / library_audit$raw_library_size
write.csv(library_audit, file.path(eda_dir, "raw_vs_tmm_effective_library_size.csv"), row.names = FALSE)

message("Drawing PCA plots")
base_pca <- function(color_var, color_label, path, continuous = FALSE) {
  plot_df <- pca_table
  p <- ggplot(plot_df, aes(x = PC1, y = PC2, label = sample_id)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_point(aes(color = .data[[color_var]], shape = sex), size = 3.4) +
    geom_text(vjust = -0.9, size = 3) +
    labs(
      title = paste0("PCA on top ", top_n, " variable genes: colored by ", color_label),
      x = sprintf("PC1 (%.1f%%)", pct_var[1]),
      y = sprintf("PC2 (%.1f%%)", pct_var[2]),
      color = color_label
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right")
  if (continuous) {
    p <- p + scale_color_gradient(low = "#2E86AB", high = "#D1495B")
  } else if (color_var == "age_group") {
    p <- p + scale_color_manual(values = c("Young" = "#2E86AB", "Old" = "#A23B72"))
  } else if (color_var == "sex") {
    p <- p + scale_color_manual(values = c("female" = "#4C78A8", "male" = "#F58518"))
  }
  ggsave(path, p, width = 8.5, height = 6.2, dpi = 180)
}

base_pca("age_group", "age group", file.path(eda_dir, "pseudobulk_pca_by_age.png"))
base_pca("sex", "sex", file.path(eda_dir, "pseudobulk_pca_by_sex.png"))
base_pca("effective_library_size", "effective library size", file.path(eda_dir, "pseudobulk_pca_by_effective_library_size.png"), continuous = TRUE)
base_pca("raw_library_size", "raw library size", file.path(eda_dir, "pseudobulk_pca_by_raw_library_size.png"), continuous = TRUE)
base_pca("cell_count", "cell count", file.path(eda_dir, "pseudobulk_pca_by_cell_count.png"), continuous = TRUE)

# Workflow-compatible primary PCA file.
file.copy(file.path(eda_dir, "pseudobulk_pca_by_age.png"), file.path(eda_dir, "pseudobulk_pca.png"), overwrite = TRUE)

message("Drawing sample correlation heatmap")
sample_cor <- cor(logcpm, method = "spearman")
write.csv(sample_cor, file.path(eda_dir, "sample_spearman_correlation.csv"))
png(file.path(eda_dir, "sample_correlation_heatmap.png"), width = 1500, height = 1300, res = 180)
heatmap(
  sample_cor,
  symm = TRUE,
  col = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
  margins = c(8, 8),
  main = "Spearman sample correlation"
)
dev.off()

message("Drawing top variable gene heatmap")
heatmap_genes <- top_variable_genes[seq_len(min(50L, length(top_variable_genes)))]
heatmap_mat <- logcpm[heatmap_genes, , drop = FALSE]
heatmap_mat_z <- t(scale(t(heatmap_mat)))
heatmap_mat_z[!is.finite(heatmap_mat_z)] <- 0
png(file.path(eda_dir, "top_variable_genes_heatmap.png"), width = 1600, height = 1600, res = 180)
heatmap(
  heatmap_mat_z,
  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  margins = c(8, 8),
  main = "Top 50 variable genes (row z-score)"
)
dev.off()

summary_table <- data.frame(
  top_variable_genes_used_for_pca = top_n,
  total_filtered_genes_available = nrow(logcpm),
  pc1_variance_percent = pct_var[1],
  pc2_variance_percent = pct_var[2],
  pc1_spearman_age_months = pc1_age_months["estimate"],
  pc1_age_group_wilcoxon_p = pc1_age_binary["p_value"],
  pc1_sex_wilcoxon_p = pc1_sex["p_value"],
  pc1_spearman_cell_count = pc1_cell_count["estimate"],
  pc1_spearman_raw_library_size = pc1_raw_library["estimate"],
  pc1_spearman_effective_library_size = pc1_effective_library["estimate"],
  extreme_pc1_and_effective_library_overlap_count = sum(pca_table$extreme_PC1_and_effective_library),
  extreme_pc1_and_raw_library_overlap_count = sum(pca_table$extreme_PC1_and_raw_library),
  stringsAsFactors = FALSE
)
write.csv(summary_table, file.path(eda_dir, "step6_eda_summary.csv"), row.names = FALSE)

report_lines <- c(
  "# Step 6: Exploratory Analysis",
  "",
  "## Inputs",
  "",
  sprintf("- TMM-normalized logCPM: `%s`", logcpm_path),
  sprintf("- Pseudobulk metadata: `%s`", metadata_path),
  sprintf("- TMM normalization factors: `%s`", norm_factors_path),
  "",
  "## What This Step Did",
  "",
  sprintf("1. Selected the top %s high-variance genes for PCA and heatmaps, rather than using all %s filtered genes.", top_n, nrow(logcpm)),
  "2. Ran PCA on the selected high-variance genes.",
  "3. Generated PCA plots colored by age group, sex, raw library size, effective library size, and cell count.",
  "4. Tested PC1 association with age, sex, cell count, raw library size, and TMM effective library size.",
  "5. Audited whether the most extreme PC1 samples also have extreme raw/effective library sizes.",
  "6. Compared raw library size with TMM effective library size.",
  "7. Generated sample-correlation and top-variable-gene heatmaps.",
  "",
  "## PCA Summary",
  "",
  sprintf("- Genes used for PCA: %s", top_n),
  sprintf("- Total filtered genes available: %s", nrow(logcpm)),
  sprintf("- PC1 variance explained: %.2f%%", pct_var[1]),
  sprintf("- PC2 variance explained: %.2f%%", pct_var[2]),
  "",
  "## PC1 Association Checks",
  "",
  sprintf("- PC1 vs age_months Spearman rho: %.3f", pc1_age_months["estimate"]),
  sprintf("- PC1 vs age_group Wilcoxon p-value: %.4g", pc1_age_binary["p_value"]),
  sprintf("- PC1 vs sex Wilcoxon p-value: %.4g", pc1_sex["p_value"]),
  sprintf("- PC1 vs cell_count Spearman rho: %.3f", pc1_cell_count["estimate"]),
  sprintf("- PC1 vs raw_library_size Spearman rho: %.3f", pc1_raw_library["estimate"]),
  sprintf("- PC1 vs effective_library_size Spearman rho: %.3f", pc1_effective_library["estimate"]),
  "",
  "## Extreme-Sample Audit",
  "",
  sprintf("- Top-2 absolute PC1 and top-2 effective-library-size overlap count: %s", sum(pca_table$extreme_PC1_and_effective_library)),
  sprintf("- Top-2 absolute PC1 and top-2 raw-library-size overlap count: %s", sum(pca_table$extreme_PC1_and_raw_library)),
  "",
  "## Sex Adjustment Note",
  "",
  "Sex is evaluated explicitly in the PCA audit. It remains partially confounded with age because both young mice are female while old mice include both female and male mice. Downstream age-effect modeling should include sex adjustment when the design matrix is estimable and should report this limitation.",
  "",
  "## Outputs",
  "",
  "- `outputs/eda/pseudobulk_pca_by_age.png`",
  "- `outputs/eda/pseudobulk_pca_by_sex.png`",
  "- `outputs/eda/pseudobulk_pca_by_effective_library_size.png`",
  "- `outputs/eda/pseudobulk_pca_by_raw_library_size.png`",
  "- `outputs/eda/pseudobulk_pca_by_cell_count.png`",
  "- `outputs/eda/sample_correlation_heatmap.png`",
  "- `outputs/eda/top_variable_genes_heatmap.png`",
  "- `outputs/eda/pc1_associations.csv`",
  "- `outputs/eda/extreme_pc1_library_audit.csv`",
  "- `outputs/eda/raw_vs_tmm_effective_library_size.csv`",
  "- `outputs/eda/top_variable_genes_for_pca.csv`"
)
writeLines(report_lines, file.path(eda_dir, "step6_eda_report.md"))

message("Step 6 complete")
message(sprintf("PCA used top %s variable genes out of %s filtered genes", top_n, nrow(logcpm)))
message(sprintf("PC1 variance: %.2f%%; PC2 variance: %.2f%%", pct_var[1], pct_var[2]))
