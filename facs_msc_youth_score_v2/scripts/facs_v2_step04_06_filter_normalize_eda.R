#!/usr/bin/env Rscript

root <- getwd()
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

out_root <- file.path(root, "outputs", "facs_v2")
processed_dir <- file.path(out_root, "processed")
qc_dir <- file.path(out_root, "qc")
eda_dir <- file.path(out_root, "eda")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(eda_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.rds")
metadata_path <- file.path(processed_dir, "facs_v2_limb_msc_mouse_metadata.csv")

cpm_threshold <- 1
min_mice_expressed <- 2
top_variable_n <- 2000

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

plot_pca <- function(scores, color_values, color_label, path, palette = NULL) {
  png(path, width = 1100, height = 850, res = 150)
  par(mar = c(5, 5, 4, 8), xpd = TRUE)
  if (is.numeric(color_values)) {
    z <- color_values
    pal <- colorRampPalette(c("#2C7BB6", "#FFFFBF", "#D7191C"))(100)
    bins <- cut(z, breaks = 100, include.lowest = TRUE, labels = FALSE)
    cols <- pal[bins]
    plot(scores$PC1, scores$PC2, pch = 19, col = cols, xlab = "PC1", ylab = "PC2", main = paste("FACS v2 PCA by", color_label))
    text(scores$PC1, scores$PC2, labels = scores$mouse, pos = 3, cex = 0.55)
    legend("right", inset = c(-0.27, 0), legend = c("low", "mid", "high"), col = pal[c(1, 50, 100)], pch = 19, bty = "n", title = color_label)
  } else {
    vals <- as.character(color_values)
    lev <- sort(unique(vals))
    if (is.null(palette)) {
      palette <- setNames(c("#B6465F", "#2878B5", "#4A7C59", "#6D597A", "#C96F3D")[seq_along(lev)], lev)
    }
    cols <- palette[vals]
    plot(scores$PC1, scores$PC2, pch = 19, col = cols, xlab = "PC1", ylab = "PC2", main = paste("FACS v2 PCA by", color_label))
    text(scores$PC1, scores$PC2, labels = scores$mouse, pos = 3, cex = 0.55)
    legend("right", inset = c(-0.27, 0), legend = names(palette), col = palette, pch = 19, bty = "n", title = color_label)
  }
  dev.off()
}

message("Reading FACS v2 pseudobulk inputs")
counts <- readRDS(counts_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(colnames(counts), metadata$mouse)) {
  stop("Pseudobulk count columns do not match mouse metadata order")
}
metadata$age_group <- factor(metadata$age_group, levels = c("Young", "Old"))
metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
metadata$age_sex_group <- interaction(metadata$sex, metadata$age_group, sep = "_", drop = TRUE)

message("Running full-data CPM gene filter for QC/EDA")
raw_library_size <- colSums(counts)
raw_cpm <- t(t(counts) / raw_library_size) * 1e6
keep <- rowSums(raw_cpm > cpm_threshold) >= min_mice_expressed
gene_qc <- data.frame(
  gene = rownames(counts),
  total_counts = rowSums(counts),
  n_mice_detected = rowSums(counts > 0),
  n_mice_cpm_gt_1 = rowSums(raw_cpm > cpm_threshold),
  keep_full_data_qc_filter = keep,
  stringsAsFactors = FALSE
)
write.csv(gene_qc, file.path(qc_dir, "step04_full_data_gene_filter_qc.csv"), row.names = FALSE)

filtered_counts <- counts[keep, , drop = FALSE]
write.csv(
  data.frame(gene = rownames(filtered_counts), filtered_counts, check.names = FALSE),
  file.path(processed_dir, "facs_v2_limb_msc_full_data_filtered_counts.csv"),
  row.names = FALSE
)
saveRDS(filtered_counts, file.path(processed_dir, "facs_v2_limb_msc_full_data_filtered_counts.rds"))

message("Calculating full-data TMM normalization for QC/EDA")
dge <- DGEList(counts = filtered_counts)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
effective_library_size <- dge$samples$lib.size * dge$samples$norm.factors
names(effective_library_size) <- rownames(dge$samples)
norm_factors <- dge$samples$norm.factors
names(norm_factors) <- rownames(dge$samples)
metadata$raw_library_size <- as.numeric(raw_library_size[metadata$mouse])
metadata$effective_library_size_tmm <- as.numeric(effective_library_size[metadata$mouse])
metadata$tmm_norm_factor <- as.numeric(norm_factors[metadata$mouse])
metadata$full_data_filtered_detected_genes <- colSums(filtered_counts > 0)
write.csv(metadata, file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_metadata_with_tmm_qc.csv"), row.names = FALSE)
write.csv(
  data.frame(gene = rownames(logcpm), logcpm, check.names = FALSE),
  file.path(processed_dir, "facs_v2_limb_msc_full_data_tmm_logcpm.csv"),
  row.names = FALSE
)
saveRDS(dge, file.path(processed_dir, "facs_v2_limb_msc_full_data_dge_tmm.rds"))
saveRDS(logcpm, file.path(processed_dir, "facs_v2_limb_msc_full_data_tmm_logcpm.rds"))

tmm_table <- data.frame(
  mouse = metadata$mouse,
  age = metadata$age,
  age_months = metadata$age_months,
  age_group = metadata$age_group,
  sex = metadata$sex,
  n_cells = metadata$n_cells,
  raw_library_size = metadata$raw_library_size,
  tmm_norm_factor = metadata$tmm_norm_factor,
  effective_library_size_tmm = metadata$effective_library_size_tmm,
  pseudobulk_detected_genes = metadata$pseudobulk_detected_genes,
  stringsAsFactors = FALSE
)
write.csv(tmm_table, file.path(qc_dir, "step05_tmm_normalization_factors.csv"), row.names = FALSE)

message("Running voom and PCA on high-variable genes")
design_factorial <- model.matrix(~ 0 + age_sex_group, data = metadata)
voom_fit <- voom(dge, design = design_factorial, plot = FALSE)
saveRDS(voom_fit, file.path(processed_dir, "facs_v2_limb_msc_full_data_voom.rds"))

gene_var <- apply(logcpm, 1, var)
top_n <- min(top_variable_n, length(gene_var))
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(top_n)]
pca <- prcomp(t(logcpm[top_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
pca_scores <- data.frame(
  mouse = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  age = metadata$age,
  age_months = metadata$age_months,
  age_group = metadata$age_group,
  sex = metadata$sex,
  age_sex_group = metadata$age_sex_group,
  n_cells = metadata$n_cells,
  raw_library_size = metadata$raw_library_size,
  effective_library_size_tmm = metadata$effective_library_size_tmm,
  pseudobulk_detected_genes = metadata$pseudobulk_detected_genes,
  stringsAsFactors = FALSE
)
write.csv(pca_scores, file.path(eda_dir, "step06_pca_scores_high_variable_genes.csv"), row.names = FALSE)
write.csv(
  data.frame(gene = top_genes, variance = gene_var[top_genes], stringsAsFactors = FALSE),
  file.path(eda_dir, "step06_top_variable_genes_for_pca.csv"),
  row.names = FALSE
)

pc_assoc <- do.call(rbind, lapply(c("PC1", "PC2", "PC3"), function(pc) {
  data.frame(
    pc = pc,
    spearman_age_months = safe_cor(pca_scores[[pc]], pca_scores$age_months),
    spearman_sex_male_numeric = safe_cor(pca_scores[[pc]], as.integer(pca_scores$sex == "male")),
    spearman_n_cells = safe_cor(pca_scores[[pc]], pca_scores$n_cells),
    spearman_raw_library_size = safe_cor(pca_scores[[pc]], pca_scores$raw_library_size),
    spearman_effective_library_size_tmm = safe_cor(pca_scores[[pc]], pca_scores$effective_library_size_tmm),
    spearman_detected_genes = safe_cor(pca_scores[[pc]], pca_scores$pseudobulk_detected_genes),
    stringsAsFactors = FALSE
  )
}))
write.csv(pc_assoc, file.path(eda_dir, "step06_pc_associations.csv"), row.names = FALSE)

plot_pca(pca_scores, pca_scores$age_group, "age group", file.path(eda_dir, "step06_pca_by_age_group.png"), c(Young = "#2C7BB6", Old = "#D7191C"))
plot_pca(pca_scores, pca_scores$sex, "sex", file.path(eda_dir, "step06_pca_by_sex.png"), c(female = "#B6465F", male = "#2878B5"))
plot_pca(pca_scores, pca_scores$age_sex_group, "age x sex", file.path(eda_dir, "step06_pca_by_age_sex_group.png"))
plot_pca(pca_scores, pca_scores$n_cells, "cell count", file.path(eda_dir, "step06_pca_by_cell_count.png"))
plot_pca(pca_scores, pca_scores$raw_library_size, "raw library size", file.path(eda_dir, "step06_pca_by_raw_library_size.png"))
plot_pca(pca_scores, pca_scores$effective_library_size_tmm, "TMM effective library size", file.path(eda_dir, "step06_pca_by_effective_library_size.png"))
plot_pca(pca_scores, pca_scores$pseudobulk_detected_genes, "detected genes", file.path(eda_dir, "step06_pca_by_detected_genes.png"))

png(file.path(qc_dir, "step05_raw_vs_tmm_effective_library_size.png"), width = 1000, height = 800, res = 150)
plot(
  metadata$raw_library_size,
  metadata$effective_library_size_tmm,
  pch = 19,
  col = ifelse(metadata$age_group == "Young", "#2C7BB6", "#D7191C"),
  xlab = "Raw pseudobulk library size",
  ylab = "TMM effective library size",
  main = "Raw vs TMM effective library size"
)
text(metadata$raw_library_size, metadata$effective_library_size_tmm, labels = metadata$mouse, pos = 3, cex = 0.55)
abline(0, 1, lty = 2, col = "#555555")
legend("topleft", legend = c("Young", "Old"), col = c("#2C7BB6", "#D7191C"), pch = 19, bty = "n")
dev.off()

message("Writing Step 04-06 report")
filter_summary <- data.frame(
  metric = c("raw_genes", "genes_retained_cpm_gt_1_in_at_least_2_mice", "genes_removed", "pca_gene_count"),
  value = c(nrow(counts), sum(keep), sum(!keep), length(top_genes))
)
write.csv(filter_summary, file.path(qc_dir, "step04_gene_filter_summary.csv"), row.names = FALSE)

report_lines <- c(
  "# FACS v2 Step 04-06: Gene Filtering, TMM Normalization, and EDA",
  "",
  "## Scope",
  "",
  "This full-data normalization is for QC and EDA only. Nested validation and model training must repeat filtering and normalization inside each training fold.",
  "",
  "No files under `data_facs` were modified.",
  "",
  "## Gene Filter",
  "",
  paste0("- Raw genes: ", nrow(counts)),
  paste0("- Retained genes: ", sum(keep)),
  paste0("- Removed genes: ", sum(!keep)),
  paste0("- Rule: CPM > ", cpm_threshold, " in at least ", min_mice_expressed, " mice"),
  "",
  "## TMM",
  "",
  paste0("- Samples: ", ncol(filtered_counts)),
  paste0("- Median raw library size: ", round(median(metadata$raw_library_size), 2)),
  paste0("- Median TMM effective library size: ", round(median(metadata$effective_library_size_tmm), 2)),
  "",
  "## PCA",
  "",
  paste0("- PCA used top variable genes: ", length(top_genes)),
  "- PCA did not use all retained genes.",
  "",
  "## PC Associations",
  "",
  paste(capture.output(print(pc_assoc, row.names = FALSE)), collapse = "\n"),
  "",
  "## Outputs",
  "",
  "- `outputs/facs_v2/processed/facs_v2_limb_msc_full_data_filtered_counts.rds`",
  "- `outputs/facs_v2/processed/facs_v2_limb_msc_full_data_tmm_logcpm.rds`",
  "- `outputs/facs_v2/processed/facs_v2_limb_msc_full_data_dge_tmm.rds`",
  "- `outputs/facs_v2/processed/facs_v2_limb_msc_full_data_voom.rds`",
  "- `outputs/facs_v2/qc/step04_full_data_gene_filter_qc.csv`",
  "- `outputs/facs_v2/qc/step05_tmm_normalization_factors.csv`",
  "- `outputs/facs_v2/eda/step06_pca_scores_high_variable_genes.csv`",
  "- `outputs/facs_v2/eda/step06_pc_associations.csv`",
  "- `outputs/facs_v2/eda/step06_pca_by_age_group.png`",
  "- `outputs/facs_v2/eda/step06_pca_by_sex.png`",
  "- `outputs/facs_v2/eda/step06_pca_by_raw_library_size.png`"
)
writeLines(report_lines, file.path(eda_dir, "step04_06_filter_normalize_eda_report.md"))
writeLines(capture.output(sessionInfo()), file.path(eda_dir, "sessionInfo_step04_06.txt"))

message("Done")
message(sprintf("Report: %s", file.path(eda_dir, "step04_06_filter_normalize_eda_report.md")))
