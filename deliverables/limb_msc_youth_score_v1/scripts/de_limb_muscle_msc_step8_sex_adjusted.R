#!/usr/bin/env Rscript

local_lib <- ".Rlibs"
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(RColorBrewer)
})

dge_path <- "data/processed/pseudobulk_dge_tmm.rds"
logcpm_path <- "data/processed/pseudobulk_logcpm.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv"
de_dir <- "outputs/de"

dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

de_out <- file.path(de_dir, "all_mice_sex_adjusted_de.csv")
sex_linked_out <- file.path(de_dir, "sex_linked_gene_audit.csv")
old_sex_out <- file.path(de_dir, "old_male_vs_old_female_effects.csv")
design_out <- file.path(de_dir, "step8_design_matrix.csv")
report_out <- file.path(de_dir, "step8_sex_adjusted_de_report.md")
volcano_out <- file.path(de_dir, "volcano_plot.png")
heatmap_out <- file.path(de_dir, "top_de_genes_heatmap.png")
mouse_plot_out <- file.path(de_dir, "top_de_genes_mouse_expression.png")

sex_linked_genes <- c("Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2")

message("Reading inputs")
dge <- readRDS(dge_path)
logcpm <- readRDS(logcpm_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(colnames(dge$counts), metadata$sample_id)) {
  stop("DGE count columns do not match labeled metadata sample_id order")
}
if (!identical(colnames(logcpm), metadata$sample_id)) {
  stop("logCPM columns do not match labeled metadata sample_id order")
}
if (!identical(rownames(dge$counts), rownames(logcpm))) {
  stop("DGE counts and logCPM rownames do not match")
}

metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
metadata$age_group <- factor(metadata$age_group, levels = c("Old", "Young"))
design <- model.matrix(~ sex + age_group, data = metadata)
if (qr(design)$rank != ncol(design)) {
  stop("Guidance-compliant design matrix is not full rank")
}
write.csv(data.frame(sample_id = metadata$sample_id, design, check.names = FALSE), design_out, row.names = FALSE)

if (any(grepl(":", colnames(design)))) {
  stop("Interaction terms are not allowed by the guidance")
}
if (any(grepl("library|cell_count", colnames(design), ignore.case = TRUE))) {
  stop("Library size/cell count must not be ordinary covariates in the DE design")
}

message("Fitting sex-adjusted edgeR quasi-likelihood model")
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design, robust = TRUE)
age_test <- glmQLFTest(fit, coef = "age_groupYoung")
sex_test <- glmQLFTest(fit, coef = "sexmale")

age_de <- topTags(age_test, n = Inf, sort.by = "none")$table
sex_de <- topTags(sex_test, n = Inf, sort.by = "none")$table

age_de$gene <- rownames(age_de)
sex_de$gene <- rownames(sex_de)

message("Adding trend and sex-sensitivity diagnostics")
age_rho <- apply(logcpm, 1, function(x) suppressWarnings(cor(x, metadata$age_months, method = "spearman")))
young_mean <- rowMeans(logcpm[, metadata$age_group == "Young", drop = FALSE])
old_mean <- rowMeans(logcpm[, metadata$age_group == "Old", drop = FALSE])

old_male <- metadata$age_group == "Old" & metadata$sex == "male"
old_female <- metadata$age_group == "Old" & metadata$sex == "female"
old_male_mean <- rowMeans(logcpm[, old_male, drop = FALSE])
old_female_mean <- rowMeans(logcpm[, old_female, drop = FALSE])
old_sex_delta <- old_male_mean - old_female_mean

old_sex_p <- apply(logcpm, 1, function(x) {
  suppressWarnings(wilcox.test(x[old_male], x[old_female], exact = FALSE)$p.value)
})

de_table <- data.frame(
  gene = age_de$gene,
  logFC = age_de$logFC,
  logCPM = age_de$logCPM,
  F = age_de$F,
  p_value = age_de$PValue,
  FDR = age_de$FDR,
  sex_adjusted_age_effect = age_de$logFC,
  direction = ifelse(age_de$logFC > 0, "young_high", "old_high"),
  mean_logcpm_young = young_mean[age_de$gene],
  mean_logcpm_old = old_mean[age_de$gene],
  spearman_logcpm_vs_age_months = age_rho[age_de$gene],
  continuous_age_trend_compatible = ifelse(
    age_de$logFC > 0,
    age_rho[age_de$gene] < 0,
    age_rho[age_de$gene] > 0
  ),
  sex_logFC_male_vs_female_adjusted_for_age = sex_de$logFC[match(age_de$gene, sex_de$gene)],
  old_male_minus_old_female_logcpm = old_sex_delta[age_de$gene],
  old_sex_wilcoxon_p = old_sex_p[age_de$gene],
  abs_old_sex_delta_over_abs_age_logFC = abs(old_sex_delta[age_de$gene]) / pmax(abs(age_de$logFC), 1e-6),
  obvious_sex_linked_gene = age_de$gene %in% sex_linked_genes,
  initial_effect_size_pass_abs_logFC_gt_0_5 = abs(age_de$logFC) > 0.5,
  initial_fdr_pass_lt_0_1 = age_de$FDR < 0.1,
  initial_candidate_for_stability_review = abs(age_de$logFC) > 0.5 &
    (ifelse(age_de$logFC > 0, age_rho[age_de$gene] < 0, age_rho[age_de$gene] > 0)) &
    !(age_de$gene %in% sex_linked_genes),
  stringsAsFactors = FALSE
)
de_table <- de_table[order(de_table$FDR, -abs(de_table$logFC), de_table$gene),]
write.csv(de_table, de_out, row.names = FALSE)

old_sex_effects <- data.frame(
  gene = rownames(logcpm),
  old_male_mean_logcpm = old_male_mean,
  old_female_mean_logcpm = old_female_mean,
  old_male_minus_old_female_logcpm = old_sex_delta,
  old_sex_wilcoxon_p = old_sex_p,
  sex_adjusted_age_logFC = de_table$logFC[match(rownames(logcpm), de_table$gene)],
  abs_old_sex_delta_over_abs_age_logFC = abs(old_sex_delta) /
    pmax(abs(de_table$logFC[match(rownames(logcpm), de_table$gene)]), 1e-6),
  obvious_sex_linked_gene = rownames(logcpm) %in% sex_linked_genes,
  stringsAsFactors = FALSE
)
old_sex_effects <- old_sex_effects[order(-abs(old_sex_effects$old_male_minus_old_female_logcpm), old_sex_effects$gene),]
write.csv(old_sex_effects, old_sex_out, row.names = FALSE)

sex_linked_audit <- data.frame(gene = sex_linked_genes, stringsAsFactors = FALSE)
sex_linked_audit$present_in_filtered_matrix <- sex_linked_audit$gene %in% rownames(logcpm)
sex_linked_audit$sex_adjusted_age_logFC <- de_table$logFC[match(sex_linked_audit$gene, de_table$gene)]
sex_linked_audit$age_FDR <- de_table$FDR[match(sex_linked_audit$gene, de_table$gene)]
sex_linked_audit$sex_logFC_male_vs_female_adjusted_for_age <- sex_de$logFC[match(sex_linked_audit$gene, sex_de$gene)]
sex_linked_audit$sex_FDR <- sex_de$FDR[match(sex_linked_audit$gene, sex_de$gene)]
sex_linked_audit$old_male_minus_old_female_logcpm <- old_sex_effects$old_male_minus_old_female_logcpm[
  match(sex_linked_audit$gene, old_sex_effects$gene)
]
sex_linked_audit$recommended_final_signature_action <- ifelse(
  sex_linked_audit$present_in_filtered_matrix,
  "exclude_from_final_signature",
  "not_present_after_filtering"
)
write.csv(sex_linked_audit, sex_linked_out, row.names = FALSE)

message("Drawing DE plots")
volcano_df <- de_table
volcano_df$neg_log10_FDR <- -log10(pmax(volcano_df$FDR, .Machine$double.xmin))
volcano_df$category <- ifelse(
  volcano_df$obvious_sex_linked_gene,
  "sex_linked",
  ifelse(volcano_df$initial_effect_size_pass_abs_logFC_gt_0_5 & volcano_df$initial_fdr_pass_lt_0_1, "initial_DE", "other")
)
p <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_FDR, color = category)) +
  geom_point(alpha = 0.65, size = 1.2) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey45", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey45", linewidth = 0.3) +
  scale_color_manual(values = c("initial_DE" = "#D1495B", "sex_linked" = "#F58518", "other" = "#777777")) +
  labs(
    title = "Sex-adjusted age-group DE: Young vs Old",
    x = "log2FC (Young / Old), adjusted for sex",
    y = "-log10(FDR)",
    color = "Category"
  ) +
  theme_classic(base_size = 12)
ggsave(volcano_out, p, width = 7.5, height = 5.8, dpi = 180)

top_genes <- head(de_table$gene[!de_table$obvious_sex_linked_gene], 30)
heatmap_mat <- logcpm[top_genes, , drop = FALSE]
heatmap_z <- t(scale(t(heatmap_mat)))
heatmap_z[!is.finite(heatmap_z)] <- 0
png(heatmap_out, width = 1500, height = 1500, res = 180)
heatmap(
  heatmap_z,
  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  margins = c(8, 9),
  main = "Top sex-adjusted age DE genes (row z-score)"
)
dev.off()

plot_genes <- head(de_table$gene[!de_table$obvious_sex_linked_gene], 12)
plot_df <- do.call(
  rbind,
  lapply(plot_genes, function(gene) {
    data.frame(
      gene = gene,
      sample_id = metadata$sample_id,
      age_months = metadata$age_months,
      age_group = metadata$age_group,
      sex = metadata$sex,
      expression = as.numeric(logcpm[gene, ]),
      stringsAsFactors = FALSE
    )
  })
)
p2 <- ggplot(plot_df, aes(x = age_months, y = expression, color = sex, shape = age_group)) +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_line(aes(group = sex), alpha = 0.25) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("female" = "#4C78A8", "male" = "#F58518")) +
  labs(
    title = "Top DE genes: mouse-level normalized expression",
    x = "Age months",
    y = "TMM log2 CPM"
  ) +
  theme_classic(base_size = 11)
ggsave(mouse_plot_out, p2, width = 12, height = 8, dpi = 180)

summary_stats <- data.frame(
  n_genes_tested = nrow(de_table),
  design_formula = "~ sex + age_group",
  design_columns = paste(colnames(design), collapse = ";"),
  design_rank = qr(design)$rank,
  design_ncol = ncol(design),
  design_full_rank = qr(design)$rank == ncol(design),
  age_contrast = "age_groupYoung (Young vs Old, adjusted for sex)",
  age_by_sex_interaction_used = FALSE,
  library_size_as_covariate_used = FALSE,
  all_12_mice_used = nrow(metadata) == 12,
  low_depth_old_female_samples_retained = all(c("18-F-50", "18-F-51") %in% metadata$sample_id),
  n_abs_logFC_gt_0_5 = sum(abs(de_table$logFC) > 0.5),
  n_FDR_lt_0_1 = sum(de_table$FDR < 0.1),
  n_abs_logFC_gt_0_5_and_FDR_lt_0_1 = sum(abs(de_table$logFC) > 0.5 & de_table$FDR < 0.1),
  n_initial_candidates_for_stability_review = sum(de_table$initial_candidate_for_stability_review),
  stringsAsFactors = FALSE
)
write.csv(summary_stats, file.path(de_dir, "step8_de_summary.csv"), row.names = FALSE)

top10 <- head(de_table, 10)
sex_linked_lines <- apply(sex_linked_audit, 1, function(row) {
  paste0("- ", row[["gene"]], ": present=", row[["present_in_filtered_matrix"]],
         ", action=", row[["recommended_final_signature_action"]])
})

report_lines <- c(
  "# Step 8: Donor-Aware Sex-Adjusted Differential Expression",
  "",
  "## Inputs",
  "",
  sprintf("- edgeR DGE object with TMM factors: `%s`", dge_path),
  sprintf("- TMM logCPM matrix: `%s`", logcpm_path),
  sprintf("- Labeled pseudobulk metadata: `%s`", metadata_path),
  "",
  "## Guidance Compliance",
  "",
  "- Primary model uses all 12 mice.",
  "- Primary model design is `~ sex + age_group`.",
  "- `Old` is the age-group reference; `age_groupYoung` is the Young-vs-Old effect.",
  "- No `sex:age_group` interaction was fitted.",
  "- Library size was not added as an ordinary covariate; edgeR uses TMM effective library size as the offset.",
  "- Low-depth old female samples `18-F-50` and `18-F-51` were retained in the primary model.",
  "",
  "## Model",
  "",
  "```r",
  "metadata$sex <- factor(metadata$sex, levels = c(\"female\", \"male\"))",
  "metadata$age_group <- factor(metadata$age_group, levels = c(\"Old\", \"Young\"))",
  "design <- model.matrix(~ sex + age_group, data = metadata)",
  "dge <- estimateDisp(dge, design)",
  "fit <- glmQLFit(dge, design, robust = TRUE)",
  "age_test <- glmQLFTest(fit, coef = \"age_groupYoung\")",
  "```",
  "",
  "## Summary",
  "",
  sprintf("- Genes tested: %s", nrow(de_table)),
  sprintf("- Design full rank: %s", qr(design)$rank == ncol(design)),
  sprintf("- Genes with |logFC| > 0.5: %s", sum(abs(de_table$logFC) > 0.5)),
  sprintf("- Genes with FDR < 0.1: %s", sum(de_table$FDR < 0.1)),
  sprintf("- Genes with |logFC| > 0.5 and FDR < 0.1: %s", sum(abs(de_table$logFC) > 0.5 & de_table$FDR < 0.1)),
  sprintf("- Initial candidates for stability review: %s", sum(de_table$initial_candidate_for_stability_review)),
  "",
  "## Sex-Linked Gene Audit",
  "",
  sex_linked_lines,
  "",
  "## Important Limitation",
  "",
  "This is a sex-adjusted model, not a sex-independent result. Young male mice are unavailable, so the model cannot test whether the young-vs-old effect behaves identically in males and females.",
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", de_out),
  sprintf("- `%s`", sex_linked_out),
  sprintf("- `%s`", old_sex_out),
  sprintf("- `%s`", volcano_out),
  sprintf("- `%s`", heatmap_out),
  sprintf("- `%s`", mouse_plot_out),
  "- `outputs/de/step8_de_summary.csv`",
  "- `outputs/de/step8_design_matrix.csv`"
)
writeLines(report_lines, report_out)

message("Step 8 complete")
message(sprintf("Genes tested: %s", nrow(de_table)))
message(sprintf("Genes |logFC| > 0.5 and FDR < 0.1: %s", sum(abs(de_table$logFC) > 0.5 & de_table$FDR < 0.1)))
