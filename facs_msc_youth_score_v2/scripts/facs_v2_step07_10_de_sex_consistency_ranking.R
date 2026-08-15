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
de_dir <- file.path(out_root, "de")
sig_dir <- file.path(out_root, "signature")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(processed_dir, "facs_v2_limb_msc_full_data_filtered_counts.rds")
metadata_path <- file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_metadata_with_tmm_qc.csv")

epsilon <- 1e-6
sex_linked_genes <- c(
  "Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2",
  "Rps4y1", "Rps4y2", "Sry", "Jarid1d"
)

add_prefix <- function(tab, prefix) {
  keep <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  tab <- tab[, intersect(keep, colnames(tab)), drop = FALSE]
  colnames(tab) <- paste0(prefix, "_", colnames(tab))
  tab
}

message("Reading full-data filtered pseudobulk inputs")
counts <- readRDS(counts_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(colnames(counts), metadata$mouse)) {
  stop("Count columns do not match metadata mouse order")
}
metadata$age_group <- factor(metadata$age_group, levels = c("Young", "Old"))
metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
metadata$age_sex_group <- factor(
  interaction(metadata$sex, metadata$age_group, sep = "_", drop = TRUE),
  levels = c("female_Young", "male_Young", "female_Old", "male_Old")
)

message("Running limma-voom age-only and factorial models")
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")

design_age_only <- model.matrix(~ age_group, data = metadata)
voom_age <- voom(dge, design_age_only, plot = FALSE)
fit_age <- eBayes(lmFit(voom_age, design_age_only))
age_tab <- topTable(fit_age, coef = "age_groupOld", number = Inf, sort.by = "none")

design_factorial <- model.matrix(~ 0 + age_sex_group, data = metadata)
voom_factorial <- voom(dge, design_factorial, plot = FALSE)
fit_factorial <- lmFit(voom_factorial, design_factorial)
contrast_matrix <- makeContrasts(
  female_age = age_sex_groupfemale_Old - age_sex_groupfemale_Young,
  male_age = age_sex_groupmale_Old - age_sex_groupmale_Young,
  common_age = ((age_sex_groupfemale_Old - age_sex_groupfemale_Young) +
    (age_sex_groupmale_Old - age_sex_groupmale_Young)) / 2,
  interaction_age_by_sex = (age_sex_groupmale_Old - age_sex_groupmale_Young) -
    (age_sex_groupfemale_Old - age_sex_groupfemale_Young),
  levels = design_factorial
)
fit_contrasts <- eBayes(contrasts.fit(fit_factorial, contrast_matrix))

female_tab <- topTable(fit_contrasts, coef = "female_age", number = Inf, sort.by = "none")
male_tab <- topTable(fit_contrasts, coef = "male_age", number = Inf, sort.by = "none")
common_tab <- topTable(fit_contrasts, coef = "common_age", number = Inf, sort.by = "none")
interaction_tab <- topTable(fit_contrasts, coef = "interaction_age_by_sex", number = Inf, sort.by = "none")

genes <- rownames(counts)
de <- data.frame(gene = genes, stringsAsFactors = FALSE)
rownames(age_tab) <- rownames(female_tab) <- rownames(male_tab) <- rownames(common_tab) <- rownames(interaction_tab) <- genes
de <- cbind(
  de,
  add_prefix(age_tab[genes, , drop = FALSE], "age_only"),
  add_prefix(female_tab[genes, , drop = FALSE], "female"),
  add_prefix(male_tab[genes, , drop = FALSE], "male"),
  add_prefix(common_tab[genes, , drop = FALSE], "common"),
  add_prefix(interaction_tab[genes, , drop = FALSE], "interaction")
)

de$sex_linked_flag <- de$gene %in% sex_linked_genes
de$female_direction <- sign(de$female_logFC)
de$male_direction <- sign(de$male_logFC)
de$common_direction <- sign(de$common_logFC)
de$age_only_direction <- sign(de$age_only_logFC)
de$sex_direction_concordant <- de$female_direction == de$male_direction & de$female_direction != 0
de$age_only_common_concordant <- de$age_only_direction == de$common_direction & de$common_direction != 0
de$interaction_magnitude_ratio <- abs(de$interaction_logFC) / (abs(de$female_logFC) + abs(de$male_logFC) + epsilon)
de$interaction_penalty <- 1 / (1 + de$interaction_magnitude_ratio)
de$primary_reliability_pass_full_data <- with(
  de,
  sex_direction_concordant &
    age_only_common_concordant &
    !sex_linked_flag &
    is.finite(common_logFC) &
    is.finite(common_t)
)
de$module <- ifelse(de$common_logFC < 0, "young_high", ifelse(de$common_logFC > 0, "old_high", "neutral"))
de$full_data_rank_score <- with(
  de,
  abs(common_logFC) * abs(common_t) * interaction_penalty *
    as.numeric(primary_reliability_pass_full_data)
)

de <- de[order(-de$full_data_rank_score, de$gene), ]
write.csv(de, file.path(de_dir, "step07_10_full_data_factorial_de_and_ranking.csv"), row.names = FALSE)

young_ranked <- de[de$module == "young_high" & de$primary_reliability_pass_full_data, ]
old_ranked <- de[de$module == "old_high" & de$primary_reliability_pass_full_data, ]
young_ranked <- young_ranked[order(-young_ranked$full_data_rank_score, young_ranked$gene), ]
old_ranked <- old_ranked[order(-old_ranked$full_data_rank_score, old_ranked$gene), ]
write.csv(young_ranked, file.path(sig_dir, "step10_full_data_young_high_ranked_candidates.csv"), row.names = FALSE)
write.csv(old_ranked, file.path(sig_dir, "step10_full_data_old_high_ranked_candidates.csv"), row.names = FALSE)
write.csv(rbind(head(young_ranked, 100), head(old_ranked, 100)), file.path(sig_dir, "step10_full_data_top100_each_direction_candidates.csv"), row.names = FALSE)

summary_table <- data.frame(
  metric = c(
    "filtered_genes_tested",
    "age_only_fdr_lt_0_1",
    "common_fdr_lt_0_1",
    "female_fdr_lt_0_1",
    "male_fdr_lt_0_1",
    "interaction_fdr_lt_0_1",
    "sex_direction_concordant_genes",
    "age_only_common_concordant_genes",
    "primary_reliability_pass_full_data",
    "young_high_ranked_candidates",
    "old_high_ranked_candidates",
    "sex_linked_genes_present"
  ),
  value = c(
    nrow(de),
    sum(de$age_only_adj.P.Val < 0.1, na.rm = TRUE),
    sum(de$common_adj.P.Val < 0.1, na.rm = TRUE),
    sum(de$female_adj.P.Val < 0.1, na.rm = TRUE),
    sum(de$male_adj.P.Val < 0.1, na.rm = TRUE),
    sum(de$interaction_adj.P.Val < 0.1, na.rm = TRUE),
    sum(de$sex_direction_concordant, na.rm = TRUE),
    sum(de$age_only_common_concordant, na.rm = TRUE),
    sum(de$primary_reliability_pass_full_data, na.rm = TRUE),
    nrow(young_ranked),
    nrow(old_ranked),
    sum(de$sex_linked_flag, na.rm = TRUE)
  )
)
write.csv(summary_table, file.path(de_dir, "step07_10_full_data_de_summary.csv"), row.names = FALSE)
write.csv(contrast_matrix, file.path(de_dir, "step07_factorial_contrast_matrix.csv"))

png(file.path(de_dir, "step07_common_vs_interaction_effects.png"), width = 1000, height = 850, res = 150)
plot(
  de$common_logFC,
  de$interaction_logFC,
  pch = 16,
  cex = 0.45,
  col = ifelse(de$primary_reliability_pass_full_data, "#2878B5", "#999999"),
  xlab = "Common age log2FC (old - young, equal-sex contrast)",
  ylab = "Age-by-sex interaction log2FC (male age effect - female age effect)",
  main = "FACS v2 full-data common age effect vs interaction"
)
abline(h = 0, v = 0, lty = 2, col = "#555555")
legend("topright", legend = c("reliability pass", "other"), col = c("#2878B5", "#999999"), pch = 16, bty = "n")
dev.off()

png(file.path(de_dir, "step07_female_vs_male_age_effects.png"), width = 1000, height = 850, res = 150)
plot(
  de$female_logFC,
  de$male_logFC,
  pch = 16,
  cex = 0.45,
  col = ifelse(de$primary_reliability_pass_full_data, "#2878B5", "#999999"),
  xlab = "Female age log2FC (old - young)",
  ylab = "Male age log2FC (old - young)",
  main = "FACS v2 full-data female vs male age effects"
)
abline(h = 0, v = 0, lty = 2, col = "#555555")
abline(0, 1, lty = 3, col = "#555555")
legend("topleft", legend = c("reliability pass", "other"), col = c("#2878B5", "#999999"), pch = 16, bty = "n")
dev.off()

report_lines <- c(
  "# FACS v2 Step 07-10: Full-Data DE, Sex Consistency, and Candidate Ranking",
  "",
  "## Scope",
  "",
  "This is a full-data exploratory analysis used to audit signal structure and generate provisional candidate rankings. It is not held-out validation evidence.",
  "",
  "Nested LOMO model training must repeat gene filtering, TMM, voom, DE, sex-consistency filtering, and ranking inside each training fold.",
  "",
  "No files under `data_facs` were modified.",
  "",
  "## Models",
  "",
  "- Age-only baseline DE: `~ age_group`.",
  "- Primary exploratory DE: `~ 0 + sex:age_group` cell-means factorial design.",
  "- Female age effect: female Old - female Young.",
  "- Male age effect: male Old - male Young.",
  "- Common age effect: average of female and male age effects.",
  "- Interaction: male age effect - female age effect.",
  "",
  "All logFC fields are log2 fold-changes with old - young orientation.",
  "",
  "## Candidate Reliability Rules",
  "",
  "- Female and male age effects must have the same nonzero direction.",
  "- Age-only and common effects must have the same nonzero direction.",
  "- Hard-coded obvious sex-linked genes are excluded.",
  paste0("- Interaction magnitude uses epsilon = ", epsilon, "."),
  "- Ranking score: `abs(common_log2FC) * abs(common moderated t) * interaction_penalty`.",
  "",
  "## Summary",
  "",
  paste(capture.output(print(summary_table, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation Notes",
  "",
  "- FDR-significant genes are not the sole definition of the final signature.",
  "- Candidate ranking uses effect size, moderated statistic, sex-direction concordance, age-only/common agreement, and interaction penalty.",
  "- Stability, low-depth sensitivity, and nested held-out performance are not assessed in this step.",
  "",
  "## Outputs",
  "",
  "- `outputs/facs_v2/de/step07_10_full_data_factorial_de_and_ranking.csv`",
  "- `outputs/facs_v2/de/step07_10_full_data_de_summary.csv`",
  "- `outputs/facs_v2/signature/step10_full_data_young_high_ranked_candidates.csv`",
  "- `outputs/facs_v2/signature/step10_full_data_old_high_ranked_candidates.csv`",
  "- `outputs/facs_v2/signature/step10_full_data_top100_each_direction_candidates.csv`",
  "- `outputs/facs_v2/de/step07_common_vs_interaction_effects.png`",
  "- `outputs/facs_v2/de/step07_female_vs_male_age_effects.png`"
)
writeLines(report_lines, file.path(de_dir, "step07_10_full_data_de_ranking_report.md"))
writeLines(capture.output(sessionInfo()), file.path(de_dir, "sessionInfo_step07_10.txt"))

message("Done")
message(sprintf("Report: %s", file.path(de_dir, "step07_10_full_data_de_ranking_report.md")))
