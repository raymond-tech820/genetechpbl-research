#!/usr/bin/env Rscript

# Step 20: donor-level bootstrap uncertainty for cross-assay Droplet sensitivity.
# Reads frozen cross-assay scores; does not retrain or modify data_droplet.

options(stringsAsFactors = FALSE)

root <- "/Users/strangenoah/Desktop/AI+X/Genetech/Youth_score"
in_dir <- file.path(root, "outputs", "facs_v2", "cross_assay_droplet")
out_dir <- file.path(in_dir, "bootstrap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260721)
n_boot <- as.integer(Sys.getenv("FACS_V2_CROSS_ASSAY_BOOT_N", "10000"))

scores_path <- file.path(in_dir, "cross_assay_droplet_scores_all_models.csv")
stopifnot(file.exists(scores_path))
scores <- read.csv(scores_path, check.names = FALSE)

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
    auc_young_vs_old = auc_rank(df$score, label_young),
    all_age_rho = safe_spearman(df$score, df$age_months),
    old_only_rho = safe_spearman(old_df$score, old_df$age_months),
    raw_library_rho = safe_spearman(df$score, df$raw_library_size),
    effective_library_rho = safe_spearman(df$score, df$effective_library_size),
    detected_genes_rho = safe_spearman(df$score, df$pseudobulk_detected_genes),
    cell_count_rho = safe_spearman(df$score, df$cell_count),
    young_minus_old_median = median(df$score[df$age_group == "Young"], na.rm = TRUE) - median(df$score[df$age_group == "Old"], na.rm = TRUE)
  )
}

summarize_boot <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values)) return(c(boot_median = NA_real_, ci_low = NA_real_, ci_high = NA_real_, valid_bootstrap_n = 0))
  c(
    boot_median = unname(median(values)),
    ci_low = unname(quantile(values, 0.025, names = FALSE, type = 7)),
    ci_high = unname(quantile(values, 0.975, names = FALSE, type = 7)),
    valid_bootstrap_n = length(values)
  )
}

models <- unique(scores$model_label)
bootstrap_rows <- list()
summary_rows <- list()
observed_rows <- list()

for (model in models) {
  df <- scores[scores$model_label == model, , drop = FALSE]
  df <- df[order(df$sample_id), , drop = FALSE]
  obs <- calc_metrics(df)
  observed_rows[[model]] <- data.frame(
    model = model,
    training_assay = unique(df$training_assay),
    test_assay = unique(df$test_assay),
    metric = names(obs),
    observed = as.numeric(obs),
    row.names = NULL
  )
  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(obs), dimnames = list(NULL, names(obs)))
  for (b in seq_len(n_boot)) {
    idx <- sample.int(nrow(df), nrow(df), replace = TRUE)
    boot_mat[b, ] <- calc_metrics(df[idx, , drop = FALSE])
  }
  for (metric in colnames(boot_mat)) {
    vals <- boot_mat[, metric]
    bootstrap_rows[[paste(model, metric, sep = "__")]] <- data.frame(
      model = model,
      metric = metric,
      bootstrap_iteration = seq_len(n_boot),
      value = vals,
      row.names = NULL
    )
    sm <- summarize_boot(vals)
    summary_rows[[paste(model, metric, sep = "__")]] <- data.frame(
      model = model,
      training_assay = unique(df$training_assay),
      test_assay = unique(df$test_assay),
      metric = metric,
      observed = unname(obs[metric]),
      boot_median = sm[["boot_median"]],
      ci_low = sm[["ci_low"]],
      ci_high = sm[["ci_high"]],
      valid_bootstrap_n = as.integer(sm[["valid_bootstrap_n"]]),
      total_bootstrap_n = n_boot,
      row.names = NULL
    )
  }
}

observed_df <- do.call(rbind, observed_rows)
bootstrap_df <- do.call(rbind, bootstrap_rows)
summary_df <- do.call(rbind, summary_rows)

write.csv(observed_df, file.path(out_dir, "cross_assay_droplet_bootstrap_observed_metrics.csv"), row.names = FALSE)
write.csv(bootstrap_df, file.path(out_dir, "cross_assay_droplet_bootstrap_samples.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "cross_assay_droplet_bootstrap_ci.csv"), row.names = FALSE)

primary_summary <- summary_df[summary_df$model %in% c("FACS_v2__factorial_medium_original", "Droplet_v1__Medium"), , drop = FALSE]
metric_order <- c("auc_young_vs_old", "all_age_rho", "old_only_rho", "raw_library_rho", "effective_library_rho", "detected_genes_rho", "cell_count_rho", "young_minus_old_median")
primary_summary$metric <- factor(primary_summary$metric, levels = metric_order)
primary_summary <- primary_summary[order(primary_summary$model, primary_summary$metric), ]
write.csv(primary_summary, file.path(out_dir, "cross_assay_droplet_primary_bootstrap_ci.csv"), row.names = FALSE)

fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
report <- c(
  "# Cross-Assay Droplet Donor Bootstrap Uncertainty",
  "",
  sprintf("Bootstrap iterations: %d donor-level resamples with replacement per model.", n_boot),
  "",
  "Purpose: quantify uncertainty for the already frozen cross-assay sensitivity results. This does not retrain models, tune thresholds, or select a new primary model.",
  "",
  "Primary focus:",
  "- `FACS_v2__factorial_medium_original`: FACS-trained frozen primary applied to Droplet.",
  "- `Droplet_v1__Medium`: Droplet-trained frozen primary applied to Droplet as same-assay reference.",
  "",
  "## Primary CI Summary",
  "",
  "| Model | Metric | Observed | Bootstrap median | 95% CI | Valid resamples |",
  "|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(primary_summary))) {
  row <- primary_summary[i, ]
  report <- c(report, sprintf(
    "| %s | %s | %s | %s | [%s, %s] | %d/%d |",
    row$model, as.character(row$metric), fmt(row$observed), fmt(row$boot_median), fmt(row$ci_low), fmt(row$ci_high), row$valid_bootstrap_n, row$total_bootstrap_n
  ))
}
report <- c(
  report,
  "",
  "## Interpretation Guardrails",
  "",
  "- The Droplet cohort has only 12 mice, including 2 Young donors; wide and boundary-touching CIs are expected.",
  "- AUC bootstrap resamples may be invalid when a resample lacks Young or Old donors; valid resample counts are reported.",
  "- Old-only rho is based on 18m/21m/24m Old mice and remains a small-sample diagnostic.",
  "- Technical-variable CIs quantify association in the Droplet application; they do not prove technical independence.",
  "",
  "## Outputs",
  "",
  "- `cross_assay_droplet_bootstrap_ci.csv`",
  "- `cross_assay_droplet_primary_bootstrap_ci.csv`",
  "- `cross_assay_droplet_bootstrap_samples.csv`",
  "- `cross_assay_droplet_bootstrap_observed_metrics.csv`"
)
writeLines(report, file.path(out_dir, "cross_assay_droplet_bootstrap_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step20_cross_assay_bootstrap.txt"))

message("Done")
message("Report: ", file.path(out_dir, "cross_assay_droplet_bootstrap_report.md"))
