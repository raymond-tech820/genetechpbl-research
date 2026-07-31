#!/usr/bin/env Rscript

# FACS Youth Score v2 Step 16: donor-level bootstrap confidence intervals.
# This script reads existing nested held-out scores and writes bootstrap summaries.
# It does not read or modify raw data_facs files.

options(stringsAsFactors = FALSE)

project_dir <- "/Users/strangenoah/Desktop/AI+X/Genetech/Youth_score"
out_dir <- file.path(project_dir, "outputs", "facs_v2", "bootstrap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260721)
n_boot <- 10000L

model_scores_file <- file.path(project_dir, "outputs", "facs_v2", "model_comparison", "model_comparison_nested_scores.csv")
elastic_scores_file <- file.path(project_dir, "outputs", "facs_v2", "baselines", "elastic_net", "facs_v2_elastic_net_heldout_score_points.csv")

stopifnot(file.exists(model_scores_file))
stopifnot(file.exists(elastic_scores_file))

model_scores <- read.csv(model_scores_file, check.names = FALSE)
elastic_scores <- read.csv(elastic_scores_file, check.names = FALSE)

wanted_models <- c(
  "factorial_medium_original",
  "factorial_stability_selected",
  "age_only_de",
  "pc1_baseline"
)

standard_scores <- model_scores[model_scores$model %in% wanted_models, c(
  "model", "mouse", "age", "age_months", "age_group", "sex", "score",
  "library_size", "detected_genes", "cell_count"
)]

elastic_standard <- data.frame(
  model = "elastic_net",
  mouse = elastic_scores$mouse,
  age = elastic_scores$age,
  age_months = elastic_scores$age_months,
  age_group = elastic_scores$age_group,
  sex = elastic_scores$sex,
  score = elastic_scores$predicted_young_probability,
  library_size = elastic_scores$pseudobulk_library_size,
  detected_genes = elastic_scores$pseudobulk_detected_genes,
  cell_count = elastic_scores$n_cells,
  check.names = FALSE
)

all_scores <- rbind(standard_scores, elastic_standard)
all_scores$age_group <- as.character(all_scores$age_group)
all_scores$age_months <- as.numeric(all_scores$age_months)
all_scores$score <- as.numeric(all_scores$score)
all_scores$library_size <- as.numeric(all_scores$library_size)
all_scores$detected_genes <- as.numeric(all_scores$detected_genes)
all_scores$cell_count <- as.numeric(all_scores$cell_count)

write.csv(all_scores, file.path(out_dir, "donor_bootstrap_input_scores.csv"), row.names = FALSE)

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

auc_rank <- function(score, label_young) {
  ok <- is.finite(score) & !is.na(label_young)
  score <- score[ok]
  label_young <- label_young[ok]
  n_pos <- sum(label_young == 1L)
  n_neg <- sum(label_young == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[label_young == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

calc_metrics <- function(df) {
  label_young <- ifelse(df$age_group == "Young", 1L, 0L)
  old_df <- df[df$age_group == "Old", , drop = FALSE]
  c(
    auc = auc_rank(df$score, label_young),
    all_age_rho = safe_spearman(df$score, df$age_months),
    old_only_rho = safe_spearman(old_df$score, old_df$age_months),
    library_rho = safe_spearman(df$score, df$library_size),
    young_old_median_difference = median(df$score[df$age_group == "Young"], na.rm = TRUE) -
      median(df$score[df$age_group == "Old"], na.rm = TRUE)
  )
}

summarize_boot <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(boot_median = NA_real_, ci_low = NA_real_, ci_high = NA_real_, valid_bootstrap_n = 0))
  }
  c(
    boot_median = unname(median(values)),
    ci_low = unname(quantile(values, 0.025, names = FALSE, type = 7)),
    ci_high = unname(quantile(values, 0.975, names = FALSE, type = 7)),
    valid_bootstrap_n = length(values)
  )
}

models <- unique(all_scores$model)
bootstrap_rows <- list()
summary_rows <- list()
observed_rows <- list()

for (m in models) {
  df <- all_scores[all_scores$model == m, , drop = FALSE]
  df <- df[order(df$mouse), , drop = FALSE]
  obs <- calc_metrics(df)
  observed_rows[[m]] <- data.frame(model = m, metric = names(obs), observed = as.numeric(obs), row.names = NULL)

  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(obs), dimnames = list(NULL, names(obs)))
  for (b in seq_len(n_boot)) {
    idx <- sample.int(nrow(df), nrow(df), replace = TRUE)
    boot_mat[b, ] <- calc_metrics(df[idx, , drop = FALSE])
  }

  for (metric in colnames(boot_mat)) {
    vals <- boot_mat[, metric]
    bootstrap_rows[[paste(m, metric, sep = "__")]] <- data.frame(
      model = m,
      metric = metric,
      bootstrap_iteration = seq_len(n_boot),
      value = vals,
      row.names = NULL
    )
    sm <- summarize_boot(vals)
    summary_rows[[paste(m, metric, sep = "__")]] <- data.frame(
      model = m,
      metric = metric,
      observed = unname(obs[metric]),
      boot_median = sm["boot_median"],
      ci_low = sm["ci_low"],
      ci_high = sm["ci_high"],
      valid_bootstrap_n = as.integer(sm["valid_bootstrap_n"]),
      total_bootstrap_n = n_boot,
      row.names = NULL
    )
  }
}

observed_df <- do.call(rbind, observed_rows)
bootstrap_df <- do.call(rbind, bootstrap_rows)
summary_df <- do.call(rbind, summary_rows)

write.csv(observed_df, file.path(out_dir, "donor_bootstrap_observed_metrics.csv"), row.names = FALSE)
write.csv(bootstrap_df, file.path(out_dir, "donor_bootstrap_samples.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "donor_bootstrap_ci.csv"), row.names = FALSE)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

report_file <- file.path(out_dir, "donor_bootstrap_report.md")
con <- file(report_file, open = "wt")
writeLines(c(
  "# FACS Youth Score v2 Step 16: Donor Bootstrap Uncertainty",
  "",
  sprintf("Bootstrap iterations: %d donor-level resamples with replacement per model.", n_boot),
  "",
  "Purpose: quantify uncertainty for predefined candidate/baseline models. This analysis is not a winner-selection or technical-independence test.",
  "",
  "Scores used:",
  "- Factorial Medium Original, Stability-Selected, Age-Only, and PC1 from `outputs/facs_v2/model_comparison/model_comparison_nested_scores.csv`.",
  "- Elastic Net from `outputs/facs_v2/baselines/elastic_net/facs_v2_elastic_net_heldout_score_points.csv`, using held-out predicted Young probability as the youth-oriented score.",
  "",
  "Metrics:",
  "- AUC: Young vs Old, higher score indicates more Young-like.",
  "- all-age rho: Spearman correlation between score and age months across all mice.",
  "- old-only rho: Spearman correlation within 18m/24m Old mice only.",
  "- library rho: Spearman correlation between score and pseudobulk library size.",
  "- young-old median difference: median(Young score) - median(Old score).",
  "",
  "Because the cohort has only 14 donors, bootstrap CIs should be interpreted as a small-sample uncertainty audit, not a substitute for external validation.",
  "",
  "## Bootstrap CI Summary",
  "",
  "| Model | Metric | Observed | Bootstrap median | 95% CI | Valid resamples |",
  "|---|---:|---:|---:|---:|---:|"
), con)

for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, ]
  writeLines(sprintf(
    "| %s | %s | %s | %s | [%s, %s] | %d/%d |",
    row$model,
    row$metric,
    fmt(row$observed),
    fmt(row$boot_median),
    fmt(row$ci_low),
    fmt(row$ci_high),
    row$valid_bootstrap_n,
    row$total_bootstrap_n
  ), con)
}

writeLines(c(
  "",
  "## Interpretation Guardrails",
  "",
  "- Wide or boundary-touching intervals are expected with 14 donors.",
  "- Old-only rho has fewer informative resamples because resampled Old donors may lack both 18m and 24m age values.",
  "- Library rho is reported as a diagnostic technical association; it is not corrected by this bootstrap.",
  "- These intervals summarize the already frozen nested-score outputs and do not rerun feature selection inside each bootstrap resample."
), con)
close(con)

cat("Wrote donor bootstrap outputs to:", out_dir, "\n")
cat("Report:", report_file, "\n")
