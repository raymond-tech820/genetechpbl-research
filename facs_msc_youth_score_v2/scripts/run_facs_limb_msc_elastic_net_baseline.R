#!/usr/bin/env Rscript

set.seed(20260721)

root <- getwd()
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(BPCells)
  library(glmnet)
})

input_dir <- file.path(root, "data_facs", "limb_muscle_msc")
matrix_dir <- file.path(input_dir, "expression_bpcells_young_old")
metadata_path <- file.path(input_dir, "facs_limb_muscle_msc_young_old_metadata.csv")
processed_dir <- file.path(root, "data_facs", "processed")
out_dir <- file.path(root, "outputs", "baselines", "elastic_net")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha_grid <- c(0.25, 0.5, 0.75, 1)
cpm_threshold <- 1
min_mice_expressed <- 2
log_pseudocount <- 1

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
    return(NA_real_)
  }
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

make_stratified_foldid <- function(y, max_folds = 3) {
  y <- as.integer(y)
  k <- min(max_folds, min(table(y)))
  if (!is.finite(k) || k < 2) {
    stop("Cannot create inner CV folds: fewer than 2 samples in one class")
  }
  foldid <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    idx <- sample(idx, length(idx))
    foldid[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  foldid
}

filter_training_genes <- function(counts_train) {
  lib <- colSums(counts_train)
  cpm <- t(t(counts_train) / lib) * 1e6
  rowSums(cpm > cpm_threshold) >= min_mice_expressed
}

make_logcpm <- function(counts) {
  lib <- colSums(counts)
  log2(t(t(counts) / lib) * 1e6 + log_pseudocount)
}

message("Reading FACS Limb Muscle MSC metadata")
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata$mouse <- metadata$mouse.id
metadata$age_months <- as.integer(sub("m$", "", metadata$age))
metadata$age_group <- factor(metadata$age_group, levels = c("Old", "Young"))
metadata$young_label <- as.integer(metadata$age_group == "Young")

message("Opening BPCells matrix")
mat <- open_matrix_dir(matrix_dir)
if (ncol(mat) != nrow(metadata)) {
  stop(sprintf("Matrix/metadata mismatch: matrix columns=%s metadata rows=%s", ncol(mat), nrow(metadata)))
}
if (!identical(as.character(colnames(mat)), as.character(metadata$index))) {
  stop("Matrix column names do not match metadata index order")
}

message("Building mouse-level pseudobulk counts")
mice <- unique(metadata$mouse)
mouse_meta <- do.call(rbind, lapply(mice, function(mouse_id) {
  idx <- which(metadata$mouse == mouse_id)
  data.frame(
    mouse = mouse_id,
    age = unique(metadata$age[idx]),
    age_months = unique(metadata$age_months[idx]),
    age_group = unique(as.character(metadata$age_group[idx])),
    sex = unique(metadata$sex[idx]),
    n_cells = length(idx),
    stringsAsFactors = FALSE
  )
}))
mouse_meta <- mouse_meta[order(mouse_meta$age_months, mouse_meta$sex, mouse_meta$mouse), ]
mice <- mouse_meta$mouse

pseudobulk <- sapply(mice, function(mouse_id) {
  idx <- which(metadata$mouse == mouse_id)
  as.numeric(rowSums(mat[, idx, drop = FALSE]))
})
rownames(pseudobulk) <- rownames(mat)
colnames(pseudobulk) <- mice
mouse_meta$pseudobulk_library_size <- as.numeric(colSums(pseudobulk[, mouse_meta$mouse, drop = FALSE]))
mouse_meta$pseudobulk_detected_genes <- as.integer(colSums(pseudobulk[, mouse_meta$mouse, drop = FALSE] > 0))
mouse_meta$young_label <- as.integer(mouse_meta$age_group == "Young")

write.csv(mouse_meta, file.path(processed_dir, "facs_limb_msc_pseudobulk_metadata.csv"), row.names = FALSE)
write.csv(
  data.frame(gene = rownames(pseudobulk), pseudobulk, check.names = FALSE),
  file.path(processed_dir, "facs_limb_msc_pseudobulk_counts.csv"),
  row.names = FALSE
)
saveRDS(pseudobulk, file.path(processed_dir, "facs_limb_msc_pseudobulk_counts.rds"))

age_sex_table <- as.data.frame.matrix(table(mouse_meta$age, mouse_meta$sex))
age_sex_table$age <- rownames(age_sex_table)
age_sex_table <- age_sex_table[, c("age", setdiff(colnames(age_sex_table), "age"))]
write.csv(age_sex_table, file.path(out_dir, "elastic_net_input_age_sex_mouse_counts.csv"), row.names = FALSE)

message("Running nested leave-one-mouse-out elastic-net baseline")
score_rows <- list()
coef_rows <- list()
tuning_rows <- list()

for (heldout_mouse in mice) {
  heldout_idx <- which(mouse_meta$mouse == heldout_mouse)
  train_idx <- setdiff(seq_len(nrow(mouse_meta)), heldout_idx)
  train_mice <- mouse_meta$mouse[train_idx]

  counts_train <- pseudobulk[, train_mice, drop = FALSE]
  counts_held <- pseudobulk[, heldout_mouse, drop = FALSE]

  keep <- filter_training_genes(counts_train)
  logcpm_train <- make_logcpm(counts_train[keep, , drop = FALSE])
  logcpm_held <- make_logcpm(counts_held[keep, , drop = FALSE])

  gene_mean <- rowMeans(logcpm_train)
  gene_sd <- apply(logcpm_train, 1, sd)
  keep_sd <- is.finite(gene_sd) & gene_sd > 0

  logcpm_train <- logcpm_train[keep_sd, , drop = FALSE]
  logcpm_held <- logcpm_held[keep_sd, , drop = FALSE]
  gene_mean <- gene_mean[keep_sd]
  gene_sd <- gene_sd[keep_sd]

  x_train <- t((logcpm_train - gene_mean) / gene_sd)
  x_held <- t((logcpm_held - gene_mean) / gene_sd)
  y_train <- mouse_meta$young_label[train_idx]

  foldid <- make_stratified_foldid(y_train, max_folds = 3)

  cv_results <- lapply(alpha_grid, function(alpha_value) {
    cvfit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      alpha = alpha_value,
      standardize = FALSE,
      foldid = foldid,
      type.measure = "deviance",
      nlambda = 100
    )
    lambda_idx <- which(cvfit$lambda == cvfit$lambda.min)[1]
    data.frame(
      alpha = alpha_value,
      lambda_min = cvfit$lambda.min,
      lambda_1se = cvfit$lambda.1se,
      cvm_min = cvfit$cvm[lambda_idx],
      cvsd_min = cvfit$cvsd[lambda_idx],
      stringsAsFactors = FALSE
    )
  })
  cv_table <- do.call(rbind, cv_results)
  cv_table$alpha_distance_from_half <- abs(cv_table$alpha - 0.5)
  best <- cv_table[order(cv_table$cvm_min, cv_table$alpha_distance_from_half), ][1, ]

  selected_lambda <- best$lambda_1se
  fit <- glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = best$alpha,
    lambda = selected_lambda,
    standardize = FALSE
  )

  pred_link <- as.numeric(predict(fit, newx = x_held, type = "link", s = selected_lambda))
  pred_prob <- as.numeric(predict(fit, newx = x_held, type = "response", s = selected_lambda))
  coef_mat <- as.matrix(coef(fit, s = selected_lambda))
  coef_df <- data.frame(
    fold_heldout_mouse = heldout_mouse,
    gene = rownames(coef_mat),
    coefficient = as.numeric(coef_mat[, 1]),
    stringsAsFactors = FALSE
  )
  coef_df <- coef_df[coef_df$coefficient != 0, , drop = FALSE]
  n_nonzero_genes <- sum(coef_df$gene != "(Intercept)")

  score_rows[[heldout_mouse]] <- data.frame(
    mouse = heldout_mouse,
    age = mouse_meta$age[heldout_idx],
    age_months = mouse_meta$age_months[heldout_idx],
    age_group = mouse_meta$age_group[heldout_idx],
    sex = mouse_meta$sex[heldout_idx],
    n_cells = mouse_meta$n_cells[heldout_idx],
    pseudobulk_library_size = mouse_meta$pseudobulk_library_size[heldout_idx],
    pseudobulk_detected_genes = mouse_meta$pseudobulk_detected_genes[heldout_idx],
    fold_type = ifelse(mouse_meta$age_group[heldout_idx] == "Young", "held_out_young", "held_out_old"),
    train_n_mice = length(train_idx),
    train_n_young = sum(mouse_meta$young_label[train_idx] == 1),
    train_n_old = sum(mouse_meta$young_label[train_idx] == 0),
    inner_cv_folds = length(unique(foldid)),
    training_gene_count_after_filter = sum(keep),
    training_gene_count_after_sd_filter = ncol(x_train),
    alpha = best$alpha,
    lambda_min = best$lambda_min,
    lambda_1se = best$lambda_1se,
    lambda_used = selected_lambda,
    inner_cv_deviance = best$cvm_min,
    inner_cv_deviance_sd = best$cvsd_min,
    n_nonzero_genes = n_nonzero_genes,
    predicted_raw_logit_score = pred_link,
    predicted_young_probability = pred_prob,
    stringsAsFactors = FALSE
  )
  coef_rows[[heldout_mouse]] <- coef_df
  tuning_rows[[heldout_mouse]] <- data.frame(
    fold_heldout_mouse = heldout_mouse,
    cv_table[, c("alpha", "lambda_min", "lambda_1se", "cvm_min", "cvsd_min")],
    selected = cv_table$alpha == best$alpha & cv_table$lambda_1se == best$lambda_1se,
    stringsAsFactors = FALSE
  )
}

scores <- do.call(rbind, score_rows)
scores <- scores[order(scores$age_months, scores$sex, scores$mouse), ]
coefs <- do.call(rbind, coef_rows)
tuning <- do.call(rbind, tuning_rows)

write.csv(scores, file.path(out_dir, "elastic_net_nested_lomo_scores.csv"), row.names = FALSE)
write.csv(coefs, file.path(out_dir, "elastic_net_nested_lomo_nonzero_coefficients.csv"), row.names = FALSE)
write.csv(tuning, file.path(out_dir, "elastic_net_nested_lomo_alpha_lambda_tuning.csv"), row.names = FALSE)

message("Summarizing baseline performance")
score <- scores$predicted_young_probability
summary_rows <- data.frame(
  metric = c(
    "n_mice",
    "n_young_mice",
    "n_old_mice",
    "spearman_all_score_vs_age_months",
    "spearman_old_only_score_vs_age_months",
    "young_median_score",
    "old_median_score",
    "young_minus_old_median_score",
    "young_old_auc",
    "spearman_score_vs_pseudobulk_library_size",
    "spearman_score_vs_pseudobulk_detected_genes",
    "spearman_score_vs_cell_count",
    "median_nonzero_genes_per_fold",
    "median_training_genes_after_filter"
  ),
  value = c(
    nrow(scores),
    sum(scores$age_group == "Young"),
    sum(scores$age_group == "Old"),
    safe_cor(score, scores$age_months),
    safe_cor(score[scores$age_group == "Old"], scores$age_months[scores$age_group == "Old"]),
    median(score[scores$age_group == "Young"]),
    median(score[scores$age_group == "Old"]),
    median(score[scores$age_group == "Young"]) - median(score[scores$age_group == "Old"]),
    auc_rank(as.integer(scores$age_group == "Young"), score),
    safe_cor(score, scores$pseudobulk_library_size),
    safe_cor(score, scores$pseudobulk_detected_genes),
    safe_cor(score, scores$n_cells),
    median(scores$n_nonzero_genes),
    median(scores$training_gene_count_after_filter)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_rows, file.path(out_dir, "elastic_net_baseline_summary.csv"), row.names = FALSE)

old_scores <- scores[scores$age_group == "Old", ]
sex_gap <- data.frame(
  slice = "old_only",
  n_female = sum(old_scores$sex == "female"),
  n_male = sum(old_scores$sex == "male"),
  female_median_score = ifelse(any(old_scores$sex == "female"), median(old_scores$predicted_young_probability[old_scores$sex == "female"]), NA_real_),
  male_median_score = ifelse(any(old_scores$sex == "male"), median(old_scores$predicted_young_probability[old_scores$sex == "male"]), NA_real_),
  male_minus_female_median = NA_real_
)
sex_gap$male_minus_female_median <- sex_gap$male_median_score - sex_gap$female_median_score
write.csv(sex_gap, file.path(out_dir, "elastic_net_old_only_sex_gap.csv"), row.names = FALSE)

png(file.path(out_dir, "elastic_net_nested_lomo_score_by_age_sex.png"), width = 1200, height = 850, res = 150)
palette_cols <- ifelse(scores$sex == "female", "#B6465F", "#2878B5")
plot(
  scores$age_months,
  scores$predicted_young_probability,
  pch = ifelse(scores$age_group == "Young", 16, 17),
  col = palette_cols,
  bg = palette_cols,
  xlab = "Age months",
  ylab = "Nested LOMO elastic-net predicted young score",
  main = "FACS Limb Muscle MSC Elastic-Net Baseline",
  ylim = c(0, 1)
)
text(scores$age_months, scores$predicted_young_probability, labels = scores$mouse, pos = 3, cex = 0.55)
legend(
  "topright",
  legend = c("female", "male", "Young", "Old"),
  col = c("#B6465F", "#2878B5", "black", "black"),
  pch = c(16, 16, 16, 17),
  bty = "n"
)
dev.off()

report_path <- file.path(out_dir, "elastic_net_baseline_report.md")
report_lines <- c(
  "# FACS Limb Muscle MSC Elastic-Net Baseline",
  "",
  "## Scope",
  "",
  "This is an internally trained baseline for the FACS Limb Muscle MSC Youth Score v2 project.",
  "",
  "It uses only the parsed TMS FACS Limb Muscle MSC pseudobulk data. It does not use GSE176206, Kaile's model, or any external labels.",
  "",
  "## Model",
  "",
  "For each outer leave-one-mouse-out fold, the held-out mouse is removed before every data-dependent step:",
  "",
  "1. training-only CPM gene filtering;",
  "2. training-only logCPM transformation;",
  "3. training-only gene mean and SD estimation;",
  "4. inner cross-validation over alpha and lambda;",
  "5. final elastic-net fit on the training mice;",
  "6. scoring of the held-out mouse.",
  "",
  "The baseline uses binomial elastic net with Young = 1 and Old = 0:",
  "",
  "\\[",
  "\\min_{\\beta_0,\\beta}",
  "\\left[-\\ell(y, \\beta_0 + X\\beta) +",
  "\\lambda\\left(\\alpha\\sum_g |\\beta_g| + \\frac{1-\\alpha}{2}\\sum_g \\beta_g^2\\right)\\right].",
  "\\]",
  "",
  sprintf("Alpha grid: `%s`; `alpha=0` ridge regression is intentionally excluded so this remains an elastic-net/LASSO baseline.", paste(alpha_grid, collapse = ", ")),
  "The final fold model uses `lambda.1se` after choosing alpha by the minimum inner-CV deviance, giving a more conservative baseline than `lambda.min`.",
  sprintf("Training-fold gene filter: CPM > %s in at least %s training mice.", cpm_threshold, min_mice_expressed),
  "",
  "The output score is `predicted_young_probability` from the nested held-out prediction. Because it is fit on only 14 mice, it should be treated as a baseline classifier score, not as a validated probability calibration.",
  "",
  "Expected glmnet warning: each outer training set contains only 13 mice, with 5-6 Young and 7-8 Old mice. Glmnet warns that one binomial class has fewer than 8 observations. This is a real small-sample limitation of the baseline, not evidence that external data or hidden labels were used.",
  "",
  "## Input audit",
  "",
  sprintf("- Mice: %s", nrow(scores)),
  sprintf("- Young mice: %s", sum(scores$age_group == "Young")),
  sprintf("- Old mice: %s", sum(scores$age_group == "Old")),
  sprintf("- Cells: %s", sum(scores$n_cells)),
  sprintf("- Genes in raw pseudobulk matrix: %s", nrow(pseudobulk)),
  "",
  "## Headline nested LOMO results",
  "",
  paste(capture.output(print(summary_rows, row.names = FALSE)), collapse = "\n"),
  "",
  "## Files",
  "",
  "- `data_facs/processed/facs_limb_msc_pseudobulk_counts.csv`",
  "- `data_facs/processed/facs_limb_msc_pseudobulk_metadata.csv`",
  "- `outputs/baselines/elastic_net/elastic_net_nested_lomo_scores.csv`",
  "- `outputs/baselines/elastic_net/elastic_net_nested_lomo_nonzero_coefficients.csv`",
  "- `outputs/baselines/elastic_net/elastic_net_nested_lomo_alpha_lambda_tuning.csv`",
  "- `outputs/baselines/elastic_net/elastic_net_baseline_summary.csv`",
  "- `outputs/baselines/elastic_net/elastic_net_nested_lomo_score_by_age_sex.png`",
  "",
  "## Not performed",
  "",
  "- GSE176206 external application: not performed because the data are not currently available.",
  "- Kaile frozen-model comparison: not performed because the frozen model file is not currently available.",
  "",
  "## Reproducibility",
  "",
  "- Seed: `20260721`",
  "- R session info: `outputs/baselines/elastic_net/sessionInfo.txt`"
)
writeLines(report_lines, report_path)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Done")
message(sprintf("Report: %s", report_path))
