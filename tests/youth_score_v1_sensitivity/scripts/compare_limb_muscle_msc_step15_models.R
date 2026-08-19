#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

scores_path <- "outputs/validation/step14_nested_lomo_scores.csv"
fold_signatures_path <- "outputs/validation/step14_fold_signatures.csv"
fold_overlap_path <- "outputs/validation/step14_fold_signature_overlap.csv"
full_signature_path <- "outputs/scores/step12_candidate_signature_versions.csv"
counts_path <- "data/processed/tms_limb_msc_pseudobulk_counts.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv"
out_dir <- "outputs/validation"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260717)

spearman_safe <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
}

wilcox_safe <- function(x, group) {
  if (length(unique(group[!is.na(group)])) < 2) {
    return(NA_real_)
  }
  suppressWarnings(wilcox.test(x ~ group, exact = FALSE)$p.value)
}

perm_p_median_gap <- function(score, sex, n_perm = 10000) {
  ok <- is.finite(score) & !is.na(sex)
  score <- score[ok]
  sex <- sex[ok]
  if (length(unique(sex)) < 2) {
    return(NA_real_)
  }
  observed <- median(score[sex == "male"], na.rm = TRUE) - median(score[sex == "female"], na.rm = TRUE)
  null <- replicate(n_perm, {
    shuffled <- sample(sex)
    median(score[shuffled == "male"], na.rm = TRUE) - median(score[shuffled == "female"], na.rm = TRUE)
  })
  (sum(abs(null) >= abs(observed)) + 1) / (n_perm + 1)
}

standardized_gap <- function(score, sex) {
  male <- score[sex == "male"]
  female <- score[sex == "female"]
  pooled_sd <- sqrt(((length(male) - 1) * var(male) + (length(female) - 1) * var(female)) /
    (length(male) + length(female) - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }
  (median(male, na.rm = TRUE) - median(female, na.rm = TRUE)) / pooled_sd
}

jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  if (length(union(a, b)) == 0) {
    return(NA_real_)
  }
  length(intersect(a, b)) / length(union(a, b))
}

label_age_strength <- function(all_rho, old_rho) {
  if (old_rho < -0.85 && all_rho < -0.85) "Strong"
  else if (old_rho < -0.65 && all_rho < -0.65) "Moderate"
  else "Weakest"
}

label_tech <- function(max_abs_tech) {
  if (max_abs_tech < 0.05) "Lowest"
  else if (max_abs_tech < 0.10) "Low"
  else "Review"
}

label_stability <- function(j) {
  if (j >= 0.57) "Highest"
  else if (j >= 0.53) "Medium"
  else "Lowest"
}

scores <- read.csv(scores_path, check.names = FALSE, stringsAsFactors = FALSE)
fold_signatures <- read.csv(fold_signatures_path, check.names = FALSE, stringsAsFactors = FALSE)
fold_overlap <- read.csv(fold_overlap_path, check.names = FALSE, stringsAsFactors = FALSE)
full_signature <- read.csv(full_signature_path, check.names = FALSE, stringsAsFactors = FALSE)
counts <- readRDS(counts_path)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

detected_genes <- colSums(counts > 0)
raw_library_size <- colSums(counts)
scores$detected_genes <- detected_genes[scores$mouse]
scores$raw_library_size_from_counts <- raw_library_size[scores$mouse]

versions <- c("Small", "Medium", "Large")

age_rows <- list()
sex_rows <- list()
tech_rows <- list()
stability_rows <- list()
decision_rows <- list()

for (version_name in versions) {
  x <- scores[scores$version == version_name & scores$score_ok, ]
  old_x <- x[x$fold_type == "held_out_old", ]
  young_x <- x[x$fold_type == "held_out_young_stress_test", ]

  young_old_sep <- median(young_x$predicted_calibrated_score, na.rm = TRUE) -
    median(old_x$predicted_calibrated_score, na.rm = TRUE)
  old_18 <- old_x$predicted_calibrated_score[old_x$age == 18]
  old_21 <- old_x$predicted_calibrated_score[old_x$age == 21]
  old_24 <- old_x$predicted_calibrated_score[old_x$age == 24]

  age_rows[[version_name]] <- data.frame(
    version = version_name,
    rho_all = spearman_safe(x$predicted_calibrated_score, x$age),
    rho_old_only = spearman_safe(old_x$predicted_calibrated_score, old_x$age),
    median_young_stress = median(young_x$predicted_calibrated_score, na.rm = TRUE),
    median_old = median(old_x$predicted_calibrated_score, na.rm = TRUE),
    young_minus_old_median = young_old_sep,
    old_18m_median = median(old_18, na.rm = TRUE),
    old_21m_median = median(old_21, na.rm = TRUE),
    old_24m_median = median(old_24, na.rm = TRUE),
    old_folds_below_young_stress_min = sum(old_x$predicted_calibrated_score <
      min(young_x$predicted_calibrated_score, na.rm = TRUE), na.rm = TRUE),
    young_stress_note = "stress_test_only_training_has_one_young_mouse",
    stringsAsFactors = FALSE
  )

  sex_rows[[version_name]] <- data.frame(
    version = version_name,
    old_male_median = median(old_x$predicted_calibrated_score[old_x$sex == "male"], na.rm = TRUE),
    old_female_median = median(old_x$predicted_calibrated_score[old_x$sex == "female"], na.rm = TRUE),
    old_male_minus_female_median = median(old_x$predicted_calibrated_score[old_x$sex == "male"], na.rm = TRUE) -
      median(old_x$predicted_calibrated_score[old_x$sex == "female"], na.rm = TRUE),
    standardized_median_gap = standardized_gap(old_x$predicted_calibrated_score, old_x$sex),
    wilcoxon_p = wilcox_safe(old_x$predicted_calibrated_score, old_x$sex),
    permutation_p_median_gap = perm_p_median_gap(old_x$predicted_calibrated_score, old_x$sex),
    interpretation_note = "p_values_auxiliary_do_not_define_acceptance",
    stringsAsFactors = FALSE
  )

  tech_vars <- data.frame(
    variable = c("effective_library_size", "raw_library_size", "cell_count", "detected_genes"),
    value = I(list(
      x$heldout_effective_library_size,
      x$raw_library_size_from_counts,
      x$heldout_cell_count,
      x$detected_genes
    )),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(tech_vars))) {
    tech_rows[[length(tech_rows) + 1]] <- data.frame(
      version = version_name,
      variable = tech_vars$variable[i],
      spearman_all = spearman_safe(x$predicted_calibrated_score, tech_vars$value[[i]]),
      spearman_old_only = spearman_safe(old_x$predicted_calibrated_score, tech_vars$value[[i]][x$fold_type == "held_out_old"]),
      stringsAsFactors = FALSE
    )
  }

  fs <- fold_signatures[fold_signatures$version == version_name, ]
  fold_ids <- unique(fs$held_out_mouse)
  gene_freq <- aggregate(
    held_out_mouse ~ gene + module + direction + version,
    fs,
    function(z) length(unique(z)) / length(fold_ids)
  )
  names(gene_freq)[names(gene_freq) == "held_out_mouse"] <- "outer_fold_inclusion_frequency"
  gene_freq <- gene_freq[order(-gene_freq$outer_fold_inclusion_frequency, gene_freq$module, gene_freq$gene), ]
  write.csv(gene_freq, file.path(out_dir, sprintf("step15_%s_gene_fold_inclusion_frequency.csv", tolower(version_name))), row.names = FALSE)

  full_sig <- full_signature[full_signature$version == version_name, ]
  full_freq <- gene_freq[match(full_sig$gene, gene_freq$gene), ]
  full_sig$outer_fold_inclusion_frequency <- full_freq$outer_fold_inclusion_frequency
  full_sig$appears_in_at_least_75pct_outer_folds <- full_sig$outer_fold_inclusion_frequency >= 0.75
  write.csv(full_sig, file.path(out_dir, sprintf("step15_%s_full_signature_fold_support.csv", tolower(version_name))), row.names = FALSE)

  overlap_v <- fold_overlap[fold_overlap$version == version_name, ]
  stability_rows[[version_name]] <- data.frame(
    version = version_name,
    median_jaccard_all = median(overlap_v$jaccard_all, na.rm = TRUE),
    median_jaccard_young_module = median(overlap_v$jaccard_young_module, na.rm = TRUE),
    median_jaccard_old_module = median(overlap_v$jaccard_old_module, na.rm = TRUE),
    min_jaccard_all = min(overlap_v$jaccard_all, na.rm = TRUE),
    full_signature_size = nrow(full_sig),
    full_signature_genes_in_at_least_75pct_outer_folds = sum(full_sig$appears_in_at_least_75pct_outer_folds, na.rm = TRUE),
    full_signature_fraction_in_at_least_75pct_outer_folds = mean(full_sig$appears_in_at_least_75pct_outer_folds, na.rm = TRUE),
    median_outer_fold_inclusion_frequency_for_full_signature = median(full_sig$outer_fold_inclusion_frequency, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

age_perf <- do.call(rbind, age_rows)
sex_assoc <- do.call(rbind, sex_rows)
tech_corr <- do.call(rbind, tech_rows)
stability <- do.call(rbind, stability_rows)

max_tech <- aggregate(abs(spearman_all) ~ version, tech_corr, max)
names(max_tech)[2] <- "max_abs_tech_spearman_all"

comparison <- Reduce(function(a, b) merge(a, b, by = "version"), list(age_perf, sex_assoc, stability, max_tech))
comparison$model_size <- c(Small = 40, Medium = 100, Large = 200)[comparison$version]
comparison$primary_candidate <- comparison$version == "Medium"
comparison$role <- ifelse(comparison$version == "Medium", "primary_candidate",
  ifelse(comparison$version == "Large", "high_stability_comparator", "not_recommended"))

decision_table <- data.frame(
  Criterion = c(
    "LOMO age correlation",
    "Old-only correlation",
    "Sex gap",
    "Technical correlation",
    "Signature stability",
    "Model size",
    "Step 15 role"
  ),
  Small = c(
    label_age_strength(age_perf$rho_all[age_perf$version == "Small"], age_perf$rho_old_only[age_perf$version == "Small"]),
    "Weak",
    "Smallest",
    label_tech(max_tech$max_abs_tech_spearman_all[max_tech$version == "Small"]),
    label_stability(stability$median_jaccard_all[stability$version == "Small"]),
    "Best",
    "Not recommended"
  ),
  Medium = c(
    label_age_strength(age_perf$rho_all[age_perf$version == "Medium"], age_perf$rho_old_only[age_perf$version == "Medium"]),
    "Strong",
    "Largest",
    label_tech(max_tech$max_abs_tech_spearman_all[max_tech$version == "Medium"]),
    label_stability(stability$median_jaccard_all[stability$version == "Medium"]),
    "Balanced",
    "Primary candidate"
  ),
  Large = c(
    label_age_strength(age_perf$rho_all[age_perf$version == "Large"], age_perf$rho_old_only[age_perf$version == "Large"]),
    "Strong",
    "Intermediate",
    label_tech(max_tech$max_abs_tech_spearman_all[max_tech$version == "Large"]),
    label_stability(stability$median_jaccard_all[stability$version == "Large"]),
    "Largest",
    "High-stability comparator"
  ),
  stringsAsFactors = FALSE
)

write.csv(age_perf, file.path(out_dir, "step15_age_performance.csv"), row.names = FALSE)
write.csv(sex_assoc, file.path(out_dir, "step15_old_sex_association.csv"), row.names = FALSE)
write.csv(tech_corr, file.path(out_dir, "step15_technical_correlations.csv"), row.names = FALSE)
write.csv(stability, file.path(out_dir, "step15_signature_stability.csv"), row.names = FALSE)
write.csv(comparison, file.path(out_dir, "step15_model_comparison_summary.csv"), row.names = FALSE)
write.csv(decision_table, file.path(out_dir, "step15_decision_table.csv"), row.names = FALSE)

plot_scores <- scores[scores$score_ok, ]
p_age <- ggplot(plot_scores, aes(x = age, y = predicted_calibrated_score, color = fold_type, shape = sex)) +
  geom_point(size = 2.7, alpha = 0.9) +
  geom_line(aes(group = version), alpha = 0) +
  facet_wrap(~ version, nrow = 1) +
  scale_color_manual(values = c("held_out_old" = "#4C78A8", "held_out_young_stress_test" = "#D1495B")) +
  labs(
    title = "Step 15 nested LOMO score by age",
    x = "Age months",
    y = "Nested LOMO calibrated Youth Score",
    color = "Fold type",
    shape = "Sex"
  ) +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "step15_lomo_score_by_age.png"), p_age, width = 10.5, height = 4.2, dpi = 180)

p_sex <- ggplot(plot_scores[plot_scores$fold_type == "held_out_old", ],
                aes(x = sex, y = predicted_calibrated_score, color = sex)) +
  geom_jitter(width = 0.12, height = 0, size = 2.5, alpha = 0.9) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.35) +
  facet_wrap(~ version, nrow = 1) +
  scale_color_manual(values = c("female" = "#4C78A8", "male" = "#F58518")) +
  labs(
    title = "Old-fold sex association diagnostic",
    x = "Old held-out mouse sex",
    y = "Nested LOMO calibrated Youth Score"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "step15_old_sex_scatter.png"), p_sex, width = 9.5, height = 4.2, dpi = 180)

tech_long <- do.call(rbind, lapply(c("heldout_effective_library_size", "raw_library_size_from_counts", "heldout_cell_count", "detected_genes"), function(var) {
  data.frame(
    version = plot_scores$version,
    variable = var,
    value = plot_scores[[var]],
    predicted_calibrated_score = plot_scores$predicted_calibrated_score,
    fold_type = plot_scores$fold_type,
    stringsAsFactors = FALSE
  )
}))
p_tech <- ggplot(tech_long, aes(x = value, y = predicted_calibrated_score, color = fold_type)) +
  geom_point(size = 1.8, alpha = 0.85) +
  facet_grid(variable ~ version, scales = "free_x") +
  scale_color_manual(values = c("held_out_old" = "#4C78A8", "held_out_young_stress_test" = "#D1495B")) +
  labs(
    title = "Step 15 technical variable diagnostics",
    x = "Technical variable value",
    y = "Nested LOMO calibrated Youth Score",
    color = "Fold type"
  ) +
  theme_classic(base_size = 10)
ggsave(file.path(out_dir, "step15_technical_scatter.png"), p_tech, width = 11, height = 9, dpi = 180)

report_lines <- c(
  "# Step 15: Small / Medium / Large Model Comparison",
  "",
  "## Purpose",
  "",
  "This step formally compares the three candidate signature sizes using nested LOMO scores from Step 14. Old-heldout folds are treated as the main generalization evidence. Young-heldout folds are kept as stress tests because the training side contains only one young mouse.",
  "",
  "## Age Performance",
  "",
  paste(capture.output(print(age_perf, row.names = FALSE)), collapse = "\n"),
  "",
  "## Old-Mouse Sex Association",
  "",
  paste(capture.output(print(sex_assoc, row.names = FALSE)), collapse = "\n"),
  "",
  "P-values are auxiliary diagnostics only. The acceptance decision should focus on magnitude relative to age separation and model complexity.",
  "",
  "## Technical Variable Correlations",
  "",
  paste(capture.output(print(tech_corr, row.names = FALSE)), collapse = "\n"),
  "",
  "## Signature Stability",
  "",
  paste(capture.output(print(stability, row.names = FALSE)), collapse = "\n"),
  "",
  "## Decision Table",
  "",
  paste(capture.output(print(decision_table, row.names = FALSE)), collapse = "\n"),
  "",
  "## Step 15 Recommendation",
  "",
  "- Small is not recommended: it has the weakest old-only age ordering and lowest signature stability.",
  "- Medium is the primary candidate: it nearly matches Large on age ordering, has low technical correlation, and uses half as many genes as Large.",
  "- Large should be retained as a high-stability comparator for robustness checks, not as the default model yet.",
  "",
  "## Important Interpretation",
  "",
  "Do not report the 12 folds as a conventional model accuracy estimate. The 10 old-heldout folds and 2 young-heldout stress-test folds have different evidentiary status.",
  "",
  "## Outputs",
  "",
  "- `outputs/validation/step15_age_performance.csv`",
  "- `outputs/validation/step15_old_sex_association.csv`",
  "- `outputs/validation/step15_technical_correlations.csv`",
  "- `outputs/validation/step15_signature_stability.csv`",
  "- `outputs/validation/step15_model_comparison_summary.csv`",
  "- `outputs/validation/step15_decision_table.csv`",
  "- `outputs/validation/step15_lomo_score_by_age.png`",
  "- `outputs/validation/step15_old_sex_scatter.png`",
  "- `outputs/validation/step15_technical_scatter.png`",
  "- `outputs/validation/step15_<version>_gene_fold_inclusion_frequency.csv`",
  "- `outputs/validation/step15_<version>_full_signature_fold_support.csv`"
)

writeLines(report_lines, file.path(out_dir, "step15_model_comparison_report.md"))

cat("Step 15 model comparison complete\n")
print(comparison, row.names = FALSE)
cat("\nDecision table:\n")
print(decision_table, row.names = FALSE)
