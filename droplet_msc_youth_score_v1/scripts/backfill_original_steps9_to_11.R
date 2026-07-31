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

counts_path <- "data/processed/pseudobulk_filtered_counts.rds"
logcpm_path <- "data/processed/pseudobulk_logcpm.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv"
full_de_path <- "outputs/de/all_mice_sex_adjusted_de.csv"
reliability_path <- "outputs/stability/gene_reliability_scores.csv"

age_trend_dir <- "outputs/age_trend"
signature_dir <- "outputs/signature"
stability_dir <- "outputs/stability"

dir.create(age_trend_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(signature_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stability_dir, recursive = TRUE, showWarnings = FALSE)

c_fc <- 0.5
sex_review_ratio_threshold <- 1
stable_old_frequency_threshold <- 0.7

message("Reading inputs")
counts <- readRDS(counts_path)
logcpm <- readRDS(logcpm_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
full_de <- read.csv(full_de_path, stringsAsFactors = FALSE, check.names = FALSE)
reliability <- read.csv(reliability_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(colnames(logcpm), metadata$sample_id)) {
  stop("logCPM columns do not match metadata")
}
if (!identical(colnames(counts), metadata$sample_id)) {
  stop("count columns do not match metadata")
}

full_de <- full_de[match(rownames(logcpm), full_de$gene),]
reliability <- reliability[match(rownames(logcpm), reliability$gene),]

fit_age_model <- function(count_matrix, sample_metadata) {
  sample_metadata$sex <- factor(sample_metadata$sex, levels = c("female", "male"))
  sample_metadata$age_group <- factor(sample_metadata$age_group, levels = c("Old", "Young"))
  design <- model.matrix(~ sex + age_group, data = sample_metadata)
  if (qr(design)$rank != ncol(design)) {
    stop("Subset design matrix is not full rank")
  }
  dge <- DGEList(counts = count_matrix, samples = sample_metadata)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "age_groupYoung")
  table <- topTags(test, n = Inf, sort.by = "none")$table
  table$gene <- rownames(table)
  subset_logcpm <- cpm(dge, log = TRUE, prior.count = 2)
  rho <- apply(subset_logcpm, 1, function(x) {
    suppressWarnings(cor(x, sample_metadata$age_months, method = "spearman"))
  })
  list(table = table, rho = rho, logcpm = subset_logcpm)
}

message("Original Step 9: continuous age trends")
spearman_p <- apply(logcpm, 1, function(x) {
  suppressWarnings(cor.test(x, metadata$age_months, method = "spearman", exact = FALSE)$p.value)
})
linear_slope <- apply(logcpm, 1, function(x) coef(lm(x ~ metadata$age_months))[2])
linear_p <- apply(logcpm, 1, function(x) summary(lm(x ~ metadata$age_months))$coefficients[2, 4])

age_trend <- data.frame(
  gene = rownames(logcpm),
  spearman_rho = full_de$spearman_logcpm_vs_age_months,
  spearman_p = spearman_p,
  linear_slope_per_month = linear_slope,
  linear_p = linear_p,
  adjusted_logFC = full_de$logFC,
  FDR = full_de$FDR,
  trend_compatible_with_binary_age_effect = full_de$continuous_age_trend_compatible,
  stringsAsFactors = FALSE
)
age_trend <- age_trend[order(-abs(age_trend$spearman_rho), age_trend$gene),]
write.csv(age_trend, file.path(age_trend_dir, "gene_age_correlations.csv"), row.names = FALSE)

plot_genes <- head(reliability$gene[reliability$passes_step9_initial_reliability], 12)
plot_df <- do.call(rbind, lapply(plot_genes, function(gene) {
  data.frame(
    gene = gene,
    sample_id = metadata$sample_id,
    age_months = metadata$age_months,
    age_group = metadata$age_group,
    sex = metadata$sex,
    expression = as.numeric(logcpm[gene,]),
    stringsAsFactors = FALSE
  )
}))
p_age <- ggplot(plot_df, aes(x = age_months, y = expression, color = sex, shape = age_group)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "grey45") +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("female" = "#4C78A8", "male" = "#F58518")) +
  labs(
    title = "Top reliable candidate trajectories by chronological age",
    x = "Age months",
    y = "TMM log2 CPM"
  ) +
  theme_classic(base_size = 11)
ggsave(file.path(age_trend_dir, "top_gene_trajectories.png"), p_age, width = 12, height = 8, dpi = 180)

message("Original Step 10: candidate aging-gene sets")
candidate_base <- data.frame(
  gene = rownames(logcpm),
  adjusted_logFC = full_de$logFC,
  FDR = full_de$FDR,
  age_rho = full_de$spearman_logcpm_vs_age_months,
  continuous_age_trend_compatible = full_de$continuous_age_trend_compatible,
  effect_size_pass_abs_logFC_gt_0_5 = abs(full_de$logFC) > c_fc,
  original_step10_candidate = abs(full_de$logFC) > c_fc & full_de$continuous_age_trend_compatible,
  direction = ifelse(full_de$logFC > 0, "young_high", "old_high"),
  LOMO_sign_rate = reliability$LOMO_sign_rate,
  no_young_mouse_reversal = !reliability$young_mouse_removal_direction_reversal,
  low_depth_sign_match = reliability$low_depth_sign_match,
  obvious_sex_linked_gene = reliability$obvious_sex_linked_gene,
  strong_sex_sensitive_gene = reliability$strong_sex_sensitive_gene,
  sex_review_flag_ratio_gt_1 = reliability$sex_effect_ratio_vs_age > sex_review_ratio_threshold,
  passes_step9_initial_reliability = reliability$passes_step9_initial_reliability,
  q_g = abs(reliability$adjusted_logFC) * abs(reliability$age_rho) * reliability$combined_reliability,
  stringsAsFactors = FALSE
)

young_candidates <- candidate_base[
  candidate_base$original_step10_candidate & candidate_base$adjusted_logFC > c_fc,
]
old_candidates <- candidate_base[
  candidate_base$original_step10_candidate & candidate_base$adjusted_logFC < -c_fc,
]
young_candidates <- young_candidates[order(-young_candidates$q_g, young_candidates$gene),]
old_candidates <- old_candidates[order(-old_candidates$q_g, old_candidates$gene),]
write.csv(young_candidates, file.path(signature_dir, "candidate_young_high_genes.csv"), row.names = FALSE)
write.csv(old_candidates, file.path(signature_dir, "candidate_old_high_genes.csv"), row.names = FALSE)

message("Original Step 11: selection frequency from LOMO candidate regeneration")
mouse_ids <- metadata$sample_id
old_mouse_ids <- metadata$sample_id[metadata$age_group == "Old"]
young_mouse_ids <- metadata$sample_id[metadata$age_group == "Young"]

selection_young <- matrix(FALSE, nrow = nrow(logcpm), ncol = length(mouse_ids), dimnames = list(rownames(logcpm), mouse_ids))
selection_old <- selection_young
subset_logfc <- matrix(NA_real_, nrow = nrow(logcpm), ncol = length(mouse_ids), dimnames = list(rownames(logcpm), mouse_ids))
subset_rho <- subset_logfc

for (mouse_id in mouse_ids) {
  message(sprintf("  regenerating candidate sets leaving out %s", mouse_id))
  keep <- metadata$sample_id != mouse_id
  fit <- fit_age_model(counts[, keep, drop = FALSE], metadata[keep, , drop = FALSE])
  genes <- fit$table$gene
  subset_logfc[genes, mouse_id] <- fit$table$logFC
  subset_rho[genes, mouse_id] <- fit$rho[genes]
  selection_young[genes, mouse_id] <- fit$table$logFC > c_fc & fit$rho[genes] < 0
  selection_old[genes, mouse_id] <- fit$table$logFC < -c_fc & fit$rho[genes] > 0
}

full_direction <- sign(full_de$logFC)
young_direction_retained <- sweep(sign(subset_logfc[, young_mouse_ids, drop = FALSE]), 1, full_direction, FUN = "==")
young_direction_retained[full_direction == 0,] <- TRUE
young_direction_retention_frequency <- rowMeans(young_direction_retained, na.rm = TRUE)
young_mouse_direction_reversal <- apply(!young_direction_retained, 1, any, na.rm = TRUE)

pi_old_young_set <- rowMeans(selection_young[, old_mouse_ids, drop = FALSE], na.rm = TRUE)
pi_old_old_set <- rowMeans(selection_old[, old_mouse_ids, drop = FALSE], na.rm = TRUE)
pi_all_young_set <- rowMeans(selection_young, na.rm = TRUE)
pi_all_old_set <- rowMeans(selection_old, na.rm = TRUE)

gene_selection_frequency <- data.frame(
  gene = rownames(logcpm),
  full_direction = ifelse(full_de$logFC > 0, "young_high", "old_high"),
  adjusted_logFC = full_de$logFC,
  age_rho = full_de$spearman_logcpm_vs_age_months,
  original_step10_candidate = candidate_base$original_step10_candidate,
  full_young_high_candidate = candidate_base$original_step10_candidate & full_de$logFC > c_fc,
  full_old_high_candidate = candidate_base$original_step10_candidate & full_de$logFC < -c_fc,
  pi_old_young_set = pi_old_young_set,
  pi_old_old_set = pi_old_old_set,
  pi_all_young_set = pi_all_young_set,
  pi_all_old_set = pi_all_old_set,
  young_direction_retention_frequency = young_direction_retention_frequency,
  young_mouse_direction_reversal = young_mouse_direction_reversal,
  LOMO_sign_rate = reliability$LOMO_sign_rate,
  max_LOMO_delta = reliability$max_LOMO_delta,
  low_depth_sign_match = reliability$low_depth_sign_match,
  obvious_sex_linked_gene = reliability$obvious_sex_linked_gene,
  strong_sex_sensitive_gene = reliability$strong_sex_sensitive_gene,
  passes_step9_initial_reliability = reliability$passes_step9_initial_reliability,
  cell_bootstrap_performed = FALSE,
  cell_bootstrap_note = "Not performed in this backfill; do not treat bootstrap stability as satisfied.",
  stringsAsFactors = FALSE
)
gene_selection_frequency$pi_old_relevant_set <- ifelse(
  gene_selection_frequency$full_young_high_candidate,
  gene_selection_frequency$pi_old_young_set,
  ifelse(gene_selection_frequency$full_old_high_candidate, gene_selection_frequency$pi_old_old_set, 0)
)
gene_selection_frequency <- gene_selection_frequency[
  order(-gene_selection_frequency$pi_old_relevant_set, !gene_selection_frequency$passes_step9_initial_reliability, -abs(gene_selection_frequency$adjusted_logFC), gene_selection_frequency$gene),
]
write.csv(gene_selection_frequency, file.path(stability_dir, "gene_selection_frequency.csv"), row.names = FALSE)

stable_base <- gene_selection_frequency[
  gene_selection_frequency$original_step10_candidate &
    gene_selection_frequency$pi_old_relevant_set >= stable_old_frequency_threshold &
    !gene_selection_frequency$young_mouse_direction_reversal &
    gene_selection_frequency$low_depth_sign_match &
    !gene_selection_frequency$obvious_sex_linked_gene &
    !gene_selection_frequency$strong_sex_sensitive_gene &
    gene_selection_frequency$passes_step9_initial_reliability,
]
stable_young <- stable_base[stable_base$full_young_high_candidate,]
stable_old <- stable_base[stable_base$full_old_high_candidate,]
stable_young <- stable_young[order(-stable_young$pi_old_relevant_set, -abs(stable_young$adjusted_logFC), stable_young$gene),]
stable_old <- stable_old[order(-stable_old$pi_old_relevant_set, -abs(stable_old$adjusted_logFC), stable_old$gene),]
write.csv(stable_young, file.path(stability_dir, "stable_young_high_genes.csv"), row.names = FALSE)
write.csv(stable_old, file.path(stability_dir, "stable_old_high_genes.csv"), row.names = FALSE)

plot_df <- rbind(
  data.frame(direction = "young_high", pi_old = gene_selection_frequency$pi_old_young_set),
  data.frame(direction = "old_high", pi_old = gene_selection_frequency$pi_old_old_set)
)
p_stab <- ggplot(plot_df, aes(x = pi_old, fill = direction)) +
  geom_histogram(binwidth = 0.1, boundary = 0, alpha = 0.75, position = "identity") +
  geom_vline(xintercept = stable_old_frequency_threshold, linetype = "dashed", color = "grey20") +
  scale_fill_manual(values = c("young_high" = "#2E86AB", "old_high" = "#A23B72")) +
  labs(
    title = "Original Step 11 LOMO old-mouse selection frequency",
    x = "Selection frequency when leaving out old mice",
    y = "Gene count"
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(stability_dir, "stability_plot.png"), p_stab, width = 7.5, height = 5.5, dpi = 180)

summary <- data.frame(
  original_step9_age_trend_genes = nrow(age_trend),
  original_step10_young_candidates = nrow(young_candidates),
  original_step10_old_candidates = nrow(old_candidates),
  original_step11_stable_young_high = nrow(stable_young),
  original_step11_stable_old_high = nrow(stable_old),
  stable_old_frequency_threshold = stable_old_frequency_threshold,
  cell_bootstrap_performed = FALSE,
  guidance_adjusted_reliability_table_used = TRUE,
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(stability_dir, "original_steps9_to_11_backfill_summary.csv"), row.names = FALSE)

report_lines <- c(
  "# Original Workflow Steps 9-11 Backfill",
  "",
  "## Purpose",
  "",
  "This backfill restores the original workflow's named intermediate products through Step 11 while preserving the later guidance constraints.",
  "",
  "## Step 9: Continuous Age Trends",
  "",
  "- Recomputed Spearman correlation between TMM logCPM and chronological age in months for every filtered gene.",
  "- Also fitted a simple linear model for expression vs age for diagnostic slope and p-value.",
  "- Output: `outputs/age_trend/gene_age_correlations.csv`.",
  "- Output: `outputs/age_trend/top_gene_trajectories.png`.",
  "",
  "## Step 10: Candidate Aging-Gene Sets",
  "",
  "- Constructed original candidate sets with `|sex-adjusted logFC| > 0.5` and compatible continuous age trend.",
  "- Young-high candidates require `logFC > 0.5` and `rho < 0`.",
  "- Old-high candidates require `logFC < -0.5` and `rho > 0`.",
  "- The files include guidance-derived columns such as LOMO stability, low-depth sign match, sex-linked flags, and sex-sensitive flags.",
  "",
  "## Step 11: Stability Selection",
  "",
  "- Regenerated candidate sets in every leave-one-mouse-out subset using `~ sex + age_group`.",
  "- Recomputed subset TMM normalization, DE logFC, and continuous-age rho for each subset.",
  "- Computed old-mouse leave-out selection frequencies for young-high and old-high sets.",
  "- Computed young-mouse direction retention as a stress test.",
  "- Stable sets require old-mouse selection frequency >= 0.7, no young-mouse direction reversal, low-depth sign match, no obvious sex-linked flag, no strong sex-sensitive flag, and Step 9 reliability pass.",
  "",
  "## What Was Not Done",
  "",
  "Cell bootstrap within donor was not performed in this backfill. The output tables explicitly mark `cell_bootstrap_performed = FALSE`; bootstrap stability should not be treated as satisfied.",
  "",
  "## Summary",
  "",
  sprintf("- Original Step 10 young-high candidates: %s", nrow(young_candidates)),
  sprintf("- Original Step 10 old-high candidates: %s", nrow(old_candidates)),
  sprintf("- Original Step 11 stable young-high genes: %s", nrow(stable_young)),
  sprintf("- Original Step 11 stable old-high genes: %s", nrow(stable_old)),
  "",
  "## Outputs",
  "",
  "- `outputs/age_trend/gene_age_correlations.csv`",
  "- `outputs/age_trend/top_gene_trajectories.png`",
  "- `outputs/signature/candidate_young_high_genes.csv`",
  "- `outputs/signature/candidate_old_high_genes.csv`",
  "- `outputs/stability/gene_selection_frequency.csv`",
  "- `outputs/stability/stable_young_high_genes.csv`",
  "- `outputs/stability/stable_old_high_genes.csv`",
  "- `outputs/stability/stability_plot.png`"
)
writeLines(report_lines, file.path(stability_dir, "original_steps9_to_11_backfill_report.md"))

message("Backfill complete")
print(summary)
