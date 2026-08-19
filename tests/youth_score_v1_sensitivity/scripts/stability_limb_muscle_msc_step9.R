#!/usr/bin/env Rscript

local_lib <- ".Rlibs"
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(edgeR)
})

counts_path <- "data/processed/pseudobulk_filtered_counts.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv"
full_de_path <- "outputs/de/all_mice_sex_adjusted_de.csv"
old_sex_path <- "outputs/de/old_male_vs_old_female_effects.csv"
sex_linked_path <- "outputs/de/sex_linked_gene_audit.csv"
stability_dir <- "outputs/stability"

dir.create(stability_dir, recursive = TRUE, showWarnings = FALSE)

lomo_out <- file.path(stability_dir, "leave_one_mouse_out_gene_stability.csv")
lomo_long_out <- file.path(stability_dir, "leave_one_mouse_out_logfc_long.csv")
low_depth_out <- file.path(stability_dir, "low_depth_exclusion_sensitivity.csv")
reliability_out <- file.path(stability_dir, "gene_reliability_scores.csv")
report_out <- file.path(stability_dir, "step9_stability_report.md")

low_depth_samples <- c("18-F-50", "18-F-51")
sex_linked_genes <- c("Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2")
sex_sensitivity_ratio_threshold <- 2
sex_sensitivity_abs_delta_threshold <- 1

message("Reading inputs")
counts <- readRDS(counts_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
full_de <- read.csv(full_de_path, stringsAsFactors = FALSE, check.names = FALSE)
old_sex <- read.csv(old_sex_path, stringsAsFactors = FALSE, check.names = FALSE)
sex_linked_audit <- read.csv(sex_linked_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(colnames(counts), metadata$sample_id)) {
  stop("Counts columns do not match metadata sample_id order")
}

fit_age_model <- function(count_matrix, sample_metadata) {
  sample_metadata$sex <- factor(sample_metadata$sex, levels = c("female", "male"))
  sample_metadata$age_group <- factor(sample_metadata$age_group, levels = c("Old", "Young"))
  design <- model.matrix(~ sex + age_group, data = sample_metadata)
  if (qr(design)$rank != ncol(design)) {
    stop(sprintf(
      "Design is not full rank for samples: %s",
      paste(sample_metadata$sample_id, collapse = ",")
    ))
  }
  dge <- DGEList(counts = count_matrix, samples = sample_metadata)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "age_groupYoung")
  table <- topTags(test, n = Inf, sort.by = "none")$table
  table$gene <- rownames(table)
  logcpm <- cpm(dge, log = TRUE, prior.count = 2)
  age_rho <- apply(logcpm, 1, function(x) {
    suppressWarnings(cor(x, sample_metadata$age_months, method = "spearman"))
  })
  list(table = table, logcpm = logcpm, age_rho = age_rho, design = design, dge = dge)
}

full_logfc <- full_de$logFC
names(full_logfc) <- full_de$gene
full_fdr <- full_de$FDR
names(full_fdr) <- full_de$gene
full_initial_candidate <- full_de$initial_candidate_for_stability_review
names(full_initial_candidate) <- full_de$gene

message("Step 9A: leave-one-mouse-out refits")
lomo_logfc <- matrix(
  NA_real_,
  nrow = nrow(counts),
  ncol = nrow(metadata),
  dimnames = list(rownames(counts), metadata$sample_id)
)
lomo_fdr <- lomo_logfc

for (sample_id in metadata$sample_id) {
  message(sprintf("  leaving out %s", sample_id))
  keep_samples <- metadata$sample_id != sample_id
  fit <- fit_age_model(counts[, keep_samples, drop = FALSE], metadata[keep_samples, , drop = FALSE])
  lomo_logfc[fit$table$gene, sample_id] <- fit$table$logFC
  lomo_fdr[fit$table$gene, sample_id] <- fit$table$FDR
}

full_sign <- sign(full_logfc[rownames(counts)])
lomo_sign <- sign(lomo_logfc)
sign_match <- sweep(lomo_sign, 1, full_sign, FUN = "==")
sign_match[full_sign == 0, ] <- TRUE
lomo_sign_rate <- rowMeans(sign_match, na.rm = TRUE)
max_lomo_delta <- apply(abs(sweep(lomo_logfc, 1, full_logfc[rownames(counts)], FUN = "-")), 1, max, na.rm = TRUE)
mean_lomo_delta <- rowMeans(abs(sweep(lomo_logfc, 1, full_logfc[rownames(counts)], FUN = "-")), na.rm = TRUE)
young_samples <- metadata$sample_id[metadata$age_group == "Young"]
old_samples <- metadata$sample_id[metadata$age_group == "Old"]
young_removal_reversal <- apply(!sign_match[, young_samples, drop = FALSE], 1, any, na.rm = TRUE)
old_removal_reversal <- apply(!sign_match[, old_samples, drop = FALSE], 1, any, na.rm = TRUE)
worst_lomo_sample <- colnames(lomo_logfc)[max.col(abs(sweep(lomo_logfc, 1, full_logfc[rownames(counts)], FUN = "-")), ties.method = "first")]
min_lomo_fdr <- apply(lomo_fdr, 1, min, na.rm = TRUE)
max_lomo_fdr <- apply(lomo_fdr, 1, max, na.rm = TRUE)

lomo_stability <- data.frame(
  gene = rownames(counts),
  adjusted_logFC = full_logfc[rownames(counts)],
  FDR = full_fdr[rownames(counts)],
  LOMO_sign_rate = lomo_sign_rate,
  max_LOMO_delta = max_lomo_delta,
  mean_LOMO_delta = mean_lomo_delta,
  worst_LOMO_sample = worst_lomo_sample,
  young_mouse_removal_direction_reversal = young_removal_reversal,
  old_mouse_removal_direction_reversal = old_removal_reversal,
  min_LOMO_FDR = min_lomo_fdr,
  max_LOMO_FDR = max_lomo_fdr,
  LOMO_pass_sign_rate_ge_0_9 = lomo_sign_rate >= 0.9,
  LOMO_pass_no_young_reversal = !young_removal_reversal,
  stringsAsFactors = FALSE
)
lomo_stability <- lomo_stability[order(!lomo_stability$LOMO_pass_sign_rate_ge_0_9, -abs(lomo_stability$adjusted_logFC), lomo_stability$gene),]
write.csv(lomo_stability, lomo_out, row.names = FALSE)

lomo_long <- do.call(
  rbind,
  lapply(colnames(lomo_logfc), function(sample_id) {
    data.frame(
      gene = rownames(lomo_logfc),
      omitted_mouse = sample_id,
      omitted_age_group = metadata$age_group[match(sample_id, metadata$sample_id)],
      omitted_sex = metadata$sex[match(sample_id, metadata$sample_id)],
      full_logFC = full_logfc[rownames(lomo_logfc)],
      LOMO_logFC = lomo_logfc[, sample_id],
      LOMO_FDR = lomo_fdr[, sample_id],
      sign_match_full = sign_match[, sample_id],
      stringsAsFactors = FALSE
    )
  })
)
write.csv(lomo_long, lomo_long_out, row.names = FALSE)

message("Step 9B: low-depth exclusion sensitivity")
low_depth_keep <- !(metadata$sample_id %in% low_depth_samples)
low_depth_fit <- fit_age_model(counts[, low_depth_keep, drop = FALSE], metadata[low_depth_keep, , drop = FALSE])
low_depth_table <- low_depth_fit$table
low_depth_age_rho <- low_depth_fit$age_rho

low_depth_logfc <- low_depth_table$logFC
names(low_depth_logfc) <- low_depth_table$gene
low_depth_fdr <- low_depth_table$FDR
names(low_depth_fdr) <- low_depth_table$gene
low_depth_trend_compatible <- ifelse(
  low_depth_logfc > 0,
  low_depth_age_rho[names(low_depth_logfc)] < 0,
  low_depth_age_rho[names(low_depth_logfc)] > 0
)
low_depth_candidate <- abs(low_depth_logfc) > 0.5 & low_depth_trend_compatible

common_genes <- rownames(counts)
low_depth_sensitivity <- data.frame(
  gene = common_genes,
  full_logFC = full_logfc[common_genes],
  reduced_logFC = low_depth_logfc[common_genes],
  full_FDR = full_fdr[common_genes],
  reduced_FDR = low_depth_fdr[common_genes],
  sign_agreement = sign(full_logfc[common_genes]) == sign(low_depth_logfc[common_genes]),
  abs_effect_size_change = abs(low_depth_logfc[common_genes] - full_logfc[common_genes]),
  full_initial_candidate = full_initial_candidate[common_genes],
  reduced_candidate_abs_logFC_gt_0_5_and_trend = low_depth_candidate[common_genes],
  candidate_status_change = full_initial_candidate[common_genes] != low_depth_candidate[common_genes],
  removed_low_depth_samples = paste(low_depth_samples, collapse = ";"),
  stringsAsFactors = FALSE
)
low_depth_sensitivity <- low_depth_sensitivity[
  order(!low_depth_sensitivity$sign_agreement, -low_depth_sensitivity$abs_effect_size_change, low_depth_sensitivity$gene),
]
write.csv(low_depth_sensitivity, low_depth_out, row.names = FALSE)
rank_correlation_full_reduced <- suppressWarnings(cor(
  full_logfc[common_genes],
  low_depth_logfc[common_genes],
  method = "spearman"
))

message("Step 9C: sex-sensitive filtering and reliability table")
old_sex_by_gene <- old_sex[match(common_genes, old_sex$gene),]
sex_linked_present <- common_genes %in% sex_linked_genes
strong_sex_sensitive <- (
  abs(old_sex_by_gene$old_male_minus_old_female_logcpm) > sex_sensitivity_abs_delta_threshold &
    old_sex_by_gene$abs_old_sex_delta_over_abs_age_logFC > sex_sensitivity_ratio_threshold
)
strong_sex_sensitive[is.na(strong_sex_sensitive)] <- FALSE

lomo_by_gene <- lomo_stability[match(common_genes, lomo_stability$gene),]
low_depth_by_gene <- low_depth_sensitivity[match(common_genes, low_depth_sensitivity$gene),]

age_rho <- full_de$spearman_logcpm_vs_age_months[match(common_genes, full_de$gene)]
trend_compatible <- full_de$continuous_age_trend_compatible[match(common_genes, full_de$gene)]

reliability <- data.frame(
  gene = common_genes,
  adjusted_logFC = full_logfc[common_genes],
  FDR = full_fdr[common_genes],
  age_rho = age_rho,
  continuous_age_trend_compatible = trend_compatible,
  LOMO_sign_rate = lomo_by_gene$LOMO_sign_rate,
  max_LOMO_delta = lomo_by_gene$max_LOMO_delta,
  young_mouse_removal_direction_reversal = lomo_by_gene$young_mouse_removal_direction_reversal,
  low_depth_sign_match = low_depth_by_gene$sign_agreement,
  low_depth_abs_effect_size_change = low_depth_by_gene$abs_effect_size_change,
  low_depth_candidate_status_change = low_depth_by_gene$candidate_status_change,
  sex_effect_old_male_minus_old_female_logcpm = old_sex_by_gene$old_male_minus_old_female_logcpm,
  sex_effect_ratio_vs_age = old_sex_by_gene$abs_old_sex_delta_over_abs_age_logFC,
  obvious_sex_linked_gene = sex_linked_present,
  strong_sex_sensitive_gene = strong_sex_sensitive,
  effect_size_pass_abs_logFC_gt_0_5 = abs(full_logfc[common_genes]) > 0.5,
  LOMO_pass = lomo_by_gene$LOMO_sign_rate >= 0.9 & !lomo_by_gene$young_mouse_removal_direction_reversal,
  low_depth_pass = low_depth_by_gene$sign_agreement,
  sex_pass = !sex_linked_present & !strong_sex_sensitive,
  exclude_reason = "",
  stringsAsFactors = FALSE
)

reliability$exclude_reason <- ifelse(!reliability$effect_size_pass_abs_logFC_gt_0_5, "weak_age_effect", reliability$exclude_reason)
reliability$exclude_reason <- ifelse(!reliability$continuous_age_trend_compatible, paste(reliability$exclude_reason, "age_trend_incompatible", sep = ";"), reliability$exclude_reason)
reliability$exclude_reason <- ifelse(!reliability$LOMO_pass, paste(reliability$exclude_reason, "LOMO_unstable", sep = ";"), reliability$exclude_reason)
reliability$exclude_reason <- ifelse(!reliability$low_depth_pass, paste(reliability$exclude_reason, "low_depth_direction_reversal", sep = ";"), reliability$exclude_reason)
reliability$exclude_reason <- ifelse(reliability$obvious_sex_linked_gene, paste(reliability$exclude_reason, "sex_linked", sep = ";"), reliability$exclude_reason)
reliability$exclude_reason <- ifelse(reliability$strong_sex_sensitive_gene, paste(reliability$exclude_reason, "sex_sensitive", sep = ";"), reliability$exclude_reason)
reliability$exclude_reason <- gsub("^;+", "", reliability$exclude_reason)

reliability$passes_step9_initial_reliability <- reliability$effect_size_pass_abs_logFC_gt_0_5 &
  reliability$continuous_age_trend_compatible &
  reliability$LOMO_pass &
  reliability$low_depth_pass &
  reliability$sex_pass

reliability$base_weight_abs_logFC <- abs(reliability$adjusted_logFC)
reliability$pi_LOMO <- pmax(0, pmin(1, reliability$LOMO_sign_rate))
reliability$pi_depth <- ifelse(reliability$low_depth_pass, 1, 0)
reliability$pi_sex <- ifelse(reliability$sex_pass, 1, 0)
reliability$combined_reliability <- reliability$pi_LOMO * reliability$pi_depth * reliability$pi_sex
reliability$stability_weighted_abs_logFC <- pmin(reliability$base_weight_abs_logFC, 3) * reliability$combined_reliability

reliability <- reliability[
  order(!reliability$passes_step9_initial_reliability, -reliability$stability_weighted_abs_logFC, reliability$gene),
]
write.csv(reliability, reliability_out, row.names = FALSE)

summary_stats <- data.frame(
  n_genes = nrow(reliability),
  low_depth_rule = "prespecified_samples_18-F-50_18-F-51",
  low_depth_samples_removed_in_sensitivity = paste(low_depth_samples, collapse = ";"),
  reduced_model_mouse_count = sum(low_depth_keep),
  reduced_model_design_full_rank = qr(low_depth_fit$design)$rank == ncol(low_depth_fit$design),
  full_reduced_logFC_spearman = rank_correlation_full_reduced,
  n_full_initial_candidates = sum(full_initial_candidate[common_genes]),
  n_LOMO_sign_rate_ge_0_9 = sum(reliability$LOMO_sign_rate >= 0.9),
  n_no_young_mouse_reversal = sum(!reliability$young_mouse_removal_direction_reversal),
  n_low_depth_sign_match = sum(reliability$low_depth_sign_match),
  n_obvious_sex_linked_present = sum(reliability$obvious_sex_linked_gene),
  n_strong_sex_sensitive = sum(reliability$strong_sex_sensitive_gene),
  n_passing_step9_initial_reliability = sum(reliability$passes_step9_initial_reliability),
  stringsAsFactors = FALSE
)
write.csv(summary_stats, file.path(stability_dir, "step9_stability_summary.csv"), row.names = FALSE)

report_lines <- c(
  "# Step 9: Stability and Robustness Filtering",
  "",
  "## Guidance Used",
  "",
  "- `outputs/de/guidance_step9.md`",
  "- `outputs/eda/guidance_review_for_next_steps.md`",
  "",
  "## What This Step Did",
  "",
  "1. Ran leave-one-mouse-out refits using the guidance-compliant model `~ sex + age_group`.",
  "2. Recorded each gene's leave-one-mouse-out logFC and FDR.",
  "3. Computed LOMO sign stability and maximum effect-size change.",
  "4. Checked whether removing either young mouse reverses the age-effect direction.",
  "5. Ran the prespecified low-depth sensitivity model after removing `18-F-50` and `18-F-51`.",
  "6. Compared full and reduced logFC values, sign agreement, effect-size change, and candidate status change.",
  "7. Used old male-versus-old female effects to mark sex-sensitive genes beyond the sex-chromosome blacklist.",
  "8. Wrote a combined gene reliability table for downstream signature construction.",
  "",
  "## Model Constraints",
  "",
  "- All refits use `~ sex + age_group` when estimable.",
  "- No age-by-sex interaction is used.",
  "- Library size is not added as a design covariate.",
  "- TMM normalization factors are recalculated inside each subset refit.",
  "",
  "## Prespecified Low-Depth Rule",
  "",
  "Based on Step 6/guidance, the low-depth sensitivity model removes the two old female low-depth samples:",
  "",
  "```text",
  "18-F-50",
  "18-F-51",
  "```",
  "",
  "## Results",
  "",
  sprintf("- Genes evaluated: %s", nrow(reliability)),
  sprintf("- Full initial candidates from Step 8: %s", sum(full_initial_candidate[common_genes])),
  sprintf("- LOMO sign-rate >= 0.9: %s", sum(reliability$LOMO_sign_rate >= 0.9)),
  sprintf("- No young-mouse removal direction reversal: %s", sum(!reliability$young_mouse_removal_direction_reversal)),
  sprintf("- Low-depth full/reduced sign match: %s", sum(reliability$low_depth_sign_match)),
  sprintf("- Strong sex-sensitive genes flagged: %s", sum(reliability$strong_sex_sensitive_gene)),
  sprintf("- Genes passing Step 9 initial reliability: %s", sum(reliability$passes_step9_initial_reliability)),
  sprintf("- Full/reduced logFC Spearman correlation: %.3f", rank_correlation_full_reduced),
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", lomo_out),
  sprintf("- `%s`", lomo_long_out),
  sprintf("- `%s`", low_depth_out),
  sprintf("- `%s`", reliability_out),
  "- `outputs/stability/step9_stability_summary.csv`"
)
writeLines(report_lines, report_out)

message("Step 9 complete")
message(sprintf("Genes passing Step 9 initial reliability: %s", sum(reliability$passes_step9_initial_reliability)))
message(sprintf("Full/reduced logFC Spearman: %.3f", rank_correlation_full_reduced))
