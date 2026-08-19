#!/usr/bin/env Rscript

# FACS Youth Score v2 Step 19: cross-assay sensitivity on TMS Droplet Limb Muscle MSC.
# Droplet data are used only for frozen scoring/evaluation, not model tuning.

options(stringsAsFactors = FALSE)

root <- "/Users/strangenoah/Desktop/AI+X/Genetech/Youth_score"
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path(root, "outputs", "facs_v2", "cross_assay_droplet")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Inputs: processed mouse-level droplet pseudobulk. Raw data_droplet is not modified.
droplet_counts_path <- file.path(root, "data_droplet", "processed", "tms_limb_msc_pseudobulk_counts.rds")
droplet_logcpm_path <- file.path(root, "data_droplet", "processed", "pseudobulk_logcpm.rds")
droplet_meta_path <- file.path(root, "data_droplet", "processed", "tms_limb_msc_pseudobulk_metadata_labeled.csv")

facs_parser_path <- file.path(root, "outputs", "facs_v2", "final_models", "parser", "score_facs_v2_youth_model.R")
facs_signature_path <- file.path(root, "outputs", "facs_v2", "final_models", "models", "facs_v2_full_data_frozen_signatures_all_models.csv")
facs_calibration_path <- file.path(root, "outputs", "facs_v2", "final_models", "models", "facs_v2_full_data_frozen_calibration.csv")

v1_parser_path <- file.path(root, "R", "score_limb_msc_youth.R")
v1_sig_medium <- file.path(root, "models", "limb_msc_general_youth_score_v1_signature.csv")
v1_sig_large <- file.path(root, "models", "limb_msc_general_youth_score_v1_large_comparator_signature.csv")
v1_sig_equal <- file.path(root, "models", "limb_msc_general_youth_score_v1_signature_equal_weight_medium.csv")
v1_calibration <- file.path(root, "models", "limb_msc_general_youth_score_v1_calibration.json")

required <- c(droplet_counts_path, droplet_logcpm_path, droplet_meta_path, facs_parser_path, facs_signature_path, facs_calibration_path, v1_parser_path, v1_sig_medium, v1_calibration)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required files:\n", paste(missing, collapse = "\n"))

droplet_counts <- readRDS(droplet_counts_path)
droplet_logcpm <- readRDS(droplet_logcpm_path)
droplet_meta <- read.csv(droplet_meta_path, check.names = FALSE)
stopifnot(identical(colnames(droplet_counts), droplet_meta$mouse_id))
stopifnot(identical(colnames(droplet_logcpm), droplet_meta$mouse_id))

if (!("cell_count" %in% names(droplet_meta))) stop("Droplet metadata lacks cell_count")
if (!("pseudobulk_detected_genes" %in% names(droplet_meta))) {
  droplet_meta$pseudobulk_detected_genes <- colSums(droplet_counts > 0)[match(droplet_meta$mouse_id, colnames(droplet_counts))]
}
if (!("raw_library_size" %in% names(droplet_meta))) {
  droplet_meta$raw_library_size <- colSums(droplet_counts)[match(droplet_meta$mouse_id, colnames(droplet_counts))]
}

source(facs_parser_path)
facs_models <- c(
  "factorial_medium_original",
  "factorial_large_original",
  "factorial_stability_selected",
  "factorial_medium_equal_weight",
  "age_only_de"
)

source(v1_parser_path)

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

auc_rank <- function(labels, scores) {
  labels <- as.integer(labels)
  ok <- is.finite(scores) & !is.na(labels)
  labels <- labels[ok]
  scores <- scores[ok]
  n_pos <- sum(labels == 1)
  n_neg <- sum(labels == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

score_summary <- function(df) {
  old <- df[df$age_group == "Old", , drop = FALSE]
  data.frame(
    model = unique(df$model_label),
    training_assay = unique(df$training_assay),
    test_assay = "Droplet",
    n_mice = nrow(df),
    auc_young_vs_old = auc_rank(as.integer(df$age_group == "Young"), df$score),
    all_age_rho = safe_cor(df$score, df$age_months),
    old_only_rho = safe_cor(old$score, old$age_months),
    young_median = median(df$score[df$age_group == "Young"], na.rm = TRUE),
    old_median = median(df$score[df$age_group == "Old"], na.rm = TRUE),
    young_minus_old_median = median(df$score[df$age_group == "Young"], na.rm = TRUE) - median(df$score[df$age_group == "Old"], na.rm = TRUE),
    raw_library_rho = safe_cor(df$score, df$raw_library_size),
    effective_library_rho = safe_cor(df$score, df$effective_library_size),
    detected_genes_rho = safe_cor(df$score, df$pseudobulk_detected_genes),
    cell_count_rho = safe_cor(df$score, df$cell_count),
    min_gene_coverage = min(df$gene_coverage, na.rm = TRUE),
    min_weighted_coverage = min(df$weighted_coverage, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

module_coverage <- function(signature, genes_available, model_label, training_assay) {
  module_col <- if ("module" %in% names(signature)) signature$module else rep("unknown", nrow(signature))
  weight <- if ("weight" %in% names(signature)) abs(signature$weight) else rep(1, nrow(signature))
  weight[!is.finite(weight)] <- 0
  rows <- lapply(split(seq_len(nrow(signature)), module_col), function(idx) {
    avail <- signature$gene[idx] %in% genes_available
    data.frame(
      model = model_label,
      training_assay = training_assay,
      module = unique(module_col[idx]),
      signature_genes = length(idx),
      covered_genes = sum(avail),
      gene_coverage = sum(avail) / length(idx),
      weighted_coverage = sum(weight[idx][avail]) / sum(weight[idx]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

missing_genes_table <- function(signature, genes_available, model_label, training_assay) {
  miss <- signature[!(signature$gene %in% genes_available), , drop = FALSE]
  if (!nrow(miss)) {
    return(data.frame(model = character(), training_assay = character(), gene = character(), module = character(), weight = numeric(), stringsAsFactors = FALSE))
  }
  data.frame(model = model_label, training_assay = training_assay, gene = miss$gene, module = miss$module, weight = if ("weight" %in% names(miss)) miss$weight else NA_real_, stringsAsFactors = FALSE)
}

message("Scoring FACS v2 frozen models on Droplet pseudobulk")
facs_scores <- lapply(facs_models, function(m) {
  z <- score_facs_v2_youth_model(droplet_counts, facs_signature_path, facs_calibration_path, model = m)
  data.frame(
    sample_id = z$sample_id,
    model_label = paste0("FACS_v2__", m),
    model_core = m,
    training_assay = "FACS_v2",
    test_assay = "Droplet",
    score = z$calibrated_score,
    raw_score = z$raw_score,
    young_module_score = z$young_module_score,
    old_module_score = z$old_module_score,
    gene_coverage = z$gene_coverage,
    weighted_coverage = z$weighted_coverage,
    stringsAsFactors = FALSE
  )
})
facs_scores <- do.call(rbind, facs_scores)

message("Scoring Droplet v1 frozen models on Droplet pseudobulk")
v1_scores <- list()
v1_scores[["Droplet_v1__Medium"]] <- score_limb_msc_youth(droplet_logcpm, signature_path = v1_sig_medium, calibration_path = v1_calibration)
if (file.exists(v1_sig_equal)) v1_scores[["Droplet_v1__Medium_equal_weight"]] <- score_limb_msc_youth(droplet_logcpm, signature_path = v1_sig_equal, calibration_path = v1_calibration, equal_weight = TRUE)

v1_scores <- do.call(rbind, lapply(names(v1_scores), function(name) {
  z <- v1_scores[[name]]
  data.frame(
    sample_id = z$sample_id,
    model_label = name,
    model_core = sub("^Droplet_v1__", "", name),
    training_assay = "Droplet_v1",
    test_assay = "Droplet",
    score = z$youth_score_raw_calibrated,
    raw_score = z$score_raw,
    young_module_score = z$score_young_module_raw,
    old_module_score = z$score_old_module_raw,
    gene_coverage = z$gene_coverage,
    weighted_coverage = z$weighted_coverage,
    stringsAsFactors = FALSE
  )
}))

all_scores <- rbind(facs_scores, v1_scores)
all_scores <- merge(all_scores, droplet_meta, by.x = "sample_id", by.y = "mouse_id", all.x = TRUE, sort = FALSE)
write.csv(all_scores, file.path(out_dir, "cross_assay_droplet_scores_all_models.csv"), row.names = FALSE)

summary_df <- do.call(rbind, lapply(split(all_scores, all_scores$model_label), score_summary))
write.csv(summary_df, file.path(out_dir, "cross_assay_droplet_model_summary.csv"), row.names = FALSE)

facs_sig <- read.csv(facs_signature_path, check.names = FALSE)
facs_cov <- do.call(rbind, lapply(facs_models, function(m) module_coverage(facs_sig[facs_sig$model == m, , drop = FALSE], rownames(droplet_counts), paste0("FACS_v2__", m), "FACS_v2")))
facs_missing <- do.call(rbind, lapply(facs_models, function(m) missing_genes_table(facs_sig[facs_sig$model == m, , drop = FALSE], rownames(droplet_counts), paste0("FACS_v2__", m), "FACS_v2")))

v1_sig_paths <- c(Droplet_v1__Medium = v1_sig_medium, Droplet_v1__Medium_equal_weight = v1_sig_equal)
v1_sig_paths <- v1_sig_paths[file.exists(v1_sig_paths)]
v1_cov <- do.call(rbind, lapply(names(v1_sig_paths), function(nm) module_coverage(read.csv(v1_sig_paths[[nm]], check.names = FALSE), rownames(droplet_logcpm), nm, "Droplet_v1")))
v1_missing <- do.call(rbind, lapply(names(v1_sig_paths), function(nm) missing_genes_table(read.csv(v1_sig_paths[[nm]], check.names = FALSE), rownames(droplet_logcpm), nm, "Droplet_v1")))
coverage_df <- rbind(facs_cov, v1_cov)
missing_df <- rbind(facs_missing, v1_missing)
write.csv(coverage_df, file.path(out_dir, "cross_assay_droplet_signature_module_coverage.csv"), row.names = FALSE)
write.csv(missing_df, file.path(out_dir, "cross_assay_droplet_missing_signature_genes.csv"), row.names = FALSE)

primary_pairs <- all_scores[all_scores$model_label %in% c("FACS_v2__factorial_medium_original", "Droplet_v1__Medium"), c("sample_id", "model_label", "score")]
primary_wide <- reshape(primary_pairs, idvar = "sample_id", timevar = "model_label", direction = "wide")
score_corr <- safe_cor(primary_wide$score.FACS_v2__factorial_medium_original, primary_wide$score.Droplet_v1__Medium)
rank_rows <- data.frame(
  sample_id = primary_wide$sample_id,
  facs_v2_medium_score = primary_wide$score.FACS_v2__factorial_medium_original,
  droplet_v1_medium_score = primary_wide$score.Droplet_v1__Medium,
  stringsAsFactors = FALSE
)
rank_rows <- merge(rank_rows, droplet_meta, by.x = "sample_id", by.y = "mouse_id", all.x = TRUE, sort = FALSE)
rank_rows$facs_v2_rank_desc <- rank(-rank_rows$facs_v2_medium_score, ties.method = "average")
rank_rows$droplet_v1_rank_desc <- rank(-rank_rows$droplet_v1_medium_score, ties.method = "average")
rank_rows$rank_difference <- rank_rows$facs_v2_rank_desc - rank_rows$droplet_v1_rank_desc
write.csv(rank_rows, file.path(out_dir, "facs_v2_vs_droplet_v1_primary_mouse_ranking.csv"), row.names = FALSE)

sig_facs_primary <- facs_sig[facs_sig$model == "factorial_medium_original", , drop = FALSE]
sig_v1_primary <- read.csv(v1_sig_medium, check.names = FALSE)
gene_overlap <- data.frame(
  comparison = "FACS_v2_factorial_medium_original_vs_Droplet_v1_Medium",
  facs_genes = nrow(sig_facs_primary),
  droplet_v1_genes = nrow(sig_v1_primary),
  overlap_genes = length(intersect(sig_facs_primary$gene, sig_v1_primary$gene)),
  union_genes = length(union(sig_facs_primary$gene, sig_v1_primary$gene)),
  jaccard = length(intersect(sig_facs_primary$gene, sig_v1_primary$gene)) / length(union(sig_facs_primary$gene, sig_v1_primary$gene)),
  score_spearman_on_droplet = score_corr,
  stringsAsFactors = FALSE
)
write.csv(gene_overlap, file.path(out_dir, "facs_v2_vs_droplet_v1_primary_overlap_and_score_correlation.csv"), row.names = FALSE)

# Plots for quick inspection.
plot_scores <- all_scores[all_scores$model_label %in% c("FACS_v2__factorial_medium_original", "Droplet_v1__Medium", "FACS_v2__factorial_large_original", "FACS_v2__factorial_stability_selected"), ]
plot_scores$model_label <- factor(plot_scores$model_label, levels = c("FACS_v2__factorial_medium_original", "FACS_v2__factorial_large_original", "FACS_v2__factorial_stability_selected", "Droplet_v1__Medium"))
p1 <- ggplot(plot_scores, aes(x = age_months, y = score, color = sex, shape = age_group)) +
  geom_hline(yintercept = c(0, 1), linewidth = 0.2, color = "grey70") +
  geom_point(size = 2.6) +
  facet_wrap(~ model_label, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(title = "Frozen youth scores on TMS Droplet limb MSC pseudobulk", x = "Age (months)", y = "Frozen calibrated score")
ggsave(file.path(fig_dir, "cross_assay_droplet_scores_by_age.png"), p1, width = 10, height = 7, dpi = 200)

p2 <- ggplot(plot_scores, aes(x = raw_library_size, y = score, color = age_group, shape = sex)) +
  geom_point(size = 2.6) +
  facet_wrap(~ model_label, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(title = "Score vs raw library size on Droplet pseudobulk", x = "Raw library size", y = "Frozen calibrated score")
ggsave(file.path(fig_dir, "cross_assay_droplet_score_vs_library_size.png"), p2, width = 10, height = 7, dpi = 200)

p3 <- ggplot(rank_rows, aes(x = droplet_v1_medium_score, y = facs_v2_medium_score, color = age_group, shape = sex)) +
  geom_point(size = 2.8) +
  theme_bw(base_size = 11) +
  labs(title = sprintf("FACS v2 vs Droplet v1 Medium scores on Droplet (Spearman %.3f)", score_corr), x = "Droplet v1 Medium score", y = "FACS v2 Medium score")
ggsave(file.path(fig_dir, "facs_v2_vs_droplet_v1_medium_score_correlation.png"), p3, width = 6, height = 5, dpi = 200)

facs_primary <- summary_df[summary_df$model == "FACS_v2__factorial_medium_original", , drop = FALSE]
droplet_primary <- summary_df[summary_df$model == "Droplet_v1__Medium", , drop = FALSE]

report <- c(
  "# FACS v2 Youth Score Cross-Assay Sensitivity on TMS Droplet Limb Muscle MSC",
  "",
  "## Scope",
  "",
  "This is a within-TMS cross-assay sensitivity / transportability analysis, not independent external validation.",
  "",
  "The FACS v2 frozen models were applied to Droplet mouse-level pseudobulk without changing signature genes, weights, training means/SDs, calibration centers, thresholds, or gene filtering. Droplet data were used only for scoring and evaluation.",
  "",
  "## Droplet Input Audit",
  "",
  sprintf("- Counts: `%s`", "data_droplet/processed/tms_limb_msc_pseudobulk_counts.rds"),
  sprintf("- Metadata: `%s`", "data_droplet/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv"),
  sprintf("- Matrix dimensions: %d genes x %d mice", nrow(droplet_counts), ncol(droplet_counts)),
  sprintf("- Mice: %d", nrow(droplet_meta)),
  sprintf("- Age composition: %s", paste(names(table(droplet_meta$age)), as.integer(table(droplet_meta$age)), sep = "=", collapse = "; ")),
  "- Sex x age_group mice:",
  paste(capture.output(print(table(droplet_meta$age_group, droplet_meta$sex))), collapse = "\n"),
  "",
  "## Signature Coverage",
  "",
  paste(capture.output(print(coverage_df, row.names = FALSE)), collapse = "\n"),
  "",
  "## Model Summary on Droplet",
  "",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"),
  "",
  "## Primary FACS v2 vs Droplet v1",
  "",
  paste(capture.output(print(gene_overlap, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  sprintf("FACS v2 primary gene coverage on Droplet was %.3f and weighted coverage was %.3f.", facs_primary$min_gene_coverage, facs_primary$min_weighted_coverage),
  sprintf("FACS v2 primary on Droplet: AUC %.3f, all-age rho %.3f, old-only rho %.3f, raw-library rho %.3f, young-minus-old median %.3f.", facs_primary$auc_young_vs_old, facs_primary$all_age_rho, facs_primary$old_only_rho, facs_primary$raw_library_rho, facs_primary$young_minus_old_median),
  sprintf("Droplet v1 Medium on Droplet: AUC %.3f, all-age rho %.3f, old-only rho %.3f, raw-library rho %.3f, young-minus-old median %.3f.", droplet_primary$auc_young_vs_old, droplet_primary$all_age_rho, droplet_primary$old_only_rho, droplet_primary$raw_library_rho, droplet_primary$young_minus_old_median),
  "",
  "Interpretation should separate cross-assay direction from technical independence. Directional transfer can be considered supportive only if coverage is adequate and young/old ordering is preserved; persistent library-size association remains a limitation.",
  "",
  "## Outputs",
  "",
  "- `cross_assay_droplet_scores_all_models.csv`",
  "- `cross_assay_droplet_model_summary.csv`",
  "- `cross_assay_droplet_signature_module_coverage.csv`",
  "- `cross_assay_droplet_missing_signature_genes.csv`",
  "- `facs_v2_vs_droplet_v1_primary_mouse_ranking.csv`",
  "- `facs_v2_vs_droplet_v1_primary_overlap_and_score_correlation.csv`",
  "- `figures/cross_assay_droplet_scores_by_age.png`",
  "- `figures/cross_assay_droplet_score_vs_library_size.png`",
  "- `figures/facs_v2_vs_droplet_v1_medium_score_correlation.png`"
)
writeLines(report, file.path(out_dir, "cross_assay_droplet_sensitivity_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step19.txt"))

message("Done")
message("Report: ", file.path(out_dir, "cross_assay_droplet_sensitivity_report.md"))
