#!/usr/bin/env Rscript

set.seed(20260721)

root <- getwd()
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

suppressPackageStartupMessages({
  library(BPCells)
  library(glmnet)
})

has_edgeR <- requireNamespace("edgeR", quietly = TRUE)

input_dir <- file.path(root, "data_facs", "limb_muscle_msc")
matrix_dir <- file.path(input_dir, "expression_bpcells_young_old")
metadata_path <- file.path(input_dir, "facs_limb_muscle_msc_young_old_metadata.csv")
out_root <- file.path(root, "outputs", "facs_v2", "baselines", "elastic_net")
processed_out <- file.path(root, "outputs", "facs_v2", "processed")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_out, recursive = TRUE, showWarnings = FALSE)

alpha_grid <- c(0.25, 0.5, 0.75, 1)
cpm_threshold <- 1
min_mice_expressed <- 2
log_pseudocount <- 1
n_boot <- 2000

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

make_stratified_foldid <- function(y, max_folds = 3) {
  y <- as.integer(y)
  k <- min(max_folds, min(table(y)))
  if (!is.finite(k) || k < 2) stop("Cannot create inner CV folds: fewer than 2 samples in one class")
  foldid <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- sample(which(y == cls))
    foldid[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  foldid
}

filter_training_genes <- function(counts_train) {
  lib <- colSums(counts_train)
  cpm <- t(t(counts_train) / lib) * 1e6
  rowSums(cpm > cpm_threshold) >= min_mice_expressed
}

make_raw_logcpm <- function(counts) {
  lib <- colSums(counts)
  log2(t(t(counts) / lib) * 1e6 + log_pseudocount)
}

make_tmm_logcpm_train <- function(counts_train) {
  dge <- edgeR::DGEList(counts = counts_train)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge, log = TRUE, prior.count = log_pseudocount, normalized.lib.sizes = TRUE)
}

make_tmm_logcpm_heldout_external <- function(counts_train, counts_held, ref_name) {
  ref_counts <- counts_train[, ref_name, drop = FALSE]
  pair <- cbind(ref_counts, counts_held)
  colnames(pair) <- c(ref_name, colnames(counts_held))
  dge_pair <- edgeR::DGEList(counts = pair)
  dge_pair <- edgeR::calcNormFactors(dge_pair, method = "TMM", refColumn = 1)
  edgeR::cpm(dge_pair, log = TRUE, prior.count = log_pseudocount, normalized.lib.sizes = TRUE)[, 2, drop = FALSE]
}

metric_summary <- function(scores) {
  score <- scores$predicted_young_probability
  data.frame(
    metric = c(
      "n_mice", "n_young_mice", "n_old_mice",
      "spearman_all_score_vs_age_months",
      "spearman_old_only_score_vs_age_months",
      "young_median_score", "old_median_score",
      "young_minus_old_median_score", "young_old_auc",
      "spearman_score_vs_pseudobulk_library_size",
      "spearman_score_vs_pseudobulk_detected_genes",
      "spearman_score_vs_cell_count",
      "median_nonzero_genes_per_fold"
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
      median(scores$n_nonzero_genes)
    ),
    stringsAsFactors = FALSE
  )
}

run_nested_en <- function(pseudobulk, mouse_meta, analysis_label, representation = c("raw_logcpm", "tmm_logcpm")) {
  representation <- match.arg(representation)
  mice <- mouse_meta$mouse
  score_rows <- list()
  coef_rows <- list()
  tuning_rows <- list()
  fold_rows <- list()
  fold_balance_rows <- list()

  for (heldout_mouse in mice) {
    heldout_idx <- which(mouse_meta$mouse == heldout_mouse)
    train_idx <- setdiff(seq_len(nrow(mouse_meta)), heldout_idx)
    train_mice <- mouse_meta$mouse[train_idx]
    counts_train <- pseudobulk[, train_mice, drop = FALSE]
    counts_held <- pseudobulk[, heldout_mouse, drop = FALSE]

    keep <- filter_training_genes(counts_train)
    counts_train <- counts_train[keep, , drop = FALSE]
    counts_held <- counts_held[keep, , drop = FALSE]

    if (representation == "raw_logcpm") {
      x_expr_train <- make_raw_logcpm(counts_train)
      x_expr_held <- make_raw_logcpm(counts_held)
      tmm_ref <- NA_character_
    } else {
      if (!has_edgeR) stop("edgeR is required for TMM-logCPM sensitivity")
      lib_train <- colSums(counts_train)
      tmm_ref <- names(which.min(abs(lib_train - median(lib_train))))
      x_expr_train <- make_tmm_logcpm_train(counts_train)
      x_expr_held <- make_tmm_logcpm_heldout_external(counts_train, counts_held, tmm_ref)
    }

    gene_mean <- rowMeans(x_expr_train)
    gene_sd <- apply(x_expr_train, 1, sd)
    keep_sd <- is.finite(gene_sd) & gene_sd > 0
    x_expr_train <- x_expr_train[keep_sd, , drop = FALSE]
    x_expr_held <- x_expr_held[keep_sd, , drop = FALSE]
    gene_mean <- gene_mean[keep_sd]
    gene_sd <- gene_sd[keep_sd]

    x_train <- t((x_expr_train - gene_mean) / gene_sd)
    x_held <- t((x_expr_held - gene_mean) / gene_sd)
    y_train <- mouse_meta$young_label[train_idx]
    foldid <- make_stratified_foldid(y_train, max_folds = 3)

    fold_rows[[heldout_mouse]] <- data.frame(
      analysis = analysis_label,
      representation = representation,
      heldout_mouse = heldout_mouse,
      train_mouse = train_mice,
      train_age_group = mouse_meta$age_group[train_idx],
      train_young_label = y_train,
      inner_fold = foldid,
      stringsAsFactors = FALSE
    )
    fold_balance_rows[[heldout_mouse]] <- do.call(rbind, lapply(sort(unique(foldid)), function(fid) {
      idx <- which(foldid == fid)
      data.frame(
        analysis = analysis_label,
        representation = representation,
        heldout_mouse = heldout_mouse,
        inner_fold = fid,
        n_train = length(idx),
        n_young = sum(y_train[idx] == 1),
        n_old = sum(y_train[idx] == 0),
        age_stratified = all(table(foldid, y_train)[, "0"] > 0) && all(table(foldid, y_train)[, "1"] > 0),
        stringsAsFactors = FALSE
      )
    }))

    cv_table <- do.call(rbind, lapply(alpha_grid, function(alpha_value) {
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
    }))
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
      analysis = analysis_label,
      representation = representation,
      fold_heldout_mouse = heldout_mouse,
      gene = rownames(coef_mat),
      coefficient = as.numeric(coef_mat[, 1]),
      stringsAsFactors = FALSE
    )
    coef_df <- coef_df[coef_df$coefficient != 0, , drop = FALSE]
    n_nonzero_genes <- sum(coef_df$gene != "(Intercept)")

    score_rows[[heldout_mouse]] <- data.frame(
      analysis = analysis_label,
      representation = representation,
      mouse = heldout_mouse,
      age = mouse_meta$age[heldout_idx],
      age_months = mouse_meta$age_months[heldout_idx],
      age_group = mouse_meta$age_group[heldout_idx],
      sex = mouse_meta$sex[heldout_idx],
      n_cells = mouse_meta$n_cells[heldout_idx],
      pseudobulk_library_size = mouse_meta$pseudobulk_library_size[heldout_idx],
      pseudobulk_detected_genes = mouse_meta$pseudobulk_detected_genes[heldout_idx],
      train_n_mice = length(train_idx),
      train_n_young = sum(y_train == 1),
      train_n_old = sum(y_train == 0),
      inner_cv_folds = length(unique(foldid)),
      inner_cv_age_stratified = all(table(foldid, y_train)[, "0"] > 0) && all(table(foldid, y_train)[, "1"] > 0),
      training_gene_count_after_filter = sum(keep),
      training_gene_count_after_sd_filter = ncol(x_train),
      alpha = best$alpha,
      lambda_min = best$lambda_min,
      lambda_1se = best$lambda_1se,
      lambda_used = selected_lambda,
      inner_cv_deviance = best$cvm_min,
      inner_cv_deviance_sd = best$cvsd_min,
      tmm_external_ref_mouse = tmm_ref,
      n_nonzero_genes = n_nonzero_genes,
      predicted_raw_logit_score = pred_link,
      predicted_young_probability = pred_prob,
      stringsAsFactors = FALSE
    )
    coef_rows[[heldout_mouse]] <- coef_df
    tuning_rows[[heldout_mouse]] <- data.frame(
      analysis = analysis_label,
      representation = representation,
      fold_heldout_mouse = heldout_mouse,
      cv_table[, c("alpha", "lambda_min", "lambda_1se", "cvm_min", "cvsd_min")],
      selected = cv_table$alpha == best$alpha & cv_table$lambda_1se == best$lambda_1se,
      stringsAsFactors = FALSE
    )
  }

  list(
    scores = do.call(rbind, score_rows),
    coefs = do.call(rbind, coef_rows),
    tuning = do.call(rbind, tuning_rows),
    inner_folds = do.call(rbind, fold_rows),
    inner_fold_balance = do.call(rbind, fold_balance_rows)
  )
}

selection_stats <- function(coefs, label) {
  genes_by_fold <- split(coefs$gene[coefs$gene != "(Intercept)"], coefs$fold_heldout_mouse[coefs$gene != "(Intercept)"])
  all_genes <- sort(unique(unlist(genes_by_fold)))
  freq <- data.frame(
    analysis = label,
    gene = all_genes,
    selected_folds = as.integer(vapply(all_genes, function(g) sum(vapply(genes_by_fold, function(x) g %in% x, logical(1))), integer(1))),
    stringsAsFactors = FALSE
  )
  freq$selection_frequency <- freq$selected_folds / length(genes_by_fold)
  freq <- freq[order(-freq$selection_frequency, freq$gene), ]

  fold_names <- names(genes_by_fold)
  jac <- list()
  k <- 1
  for (i in seq_along(fold_names)) {
    for (j in seq_along(fold_names)) {
      if (j <= i) next
      a <- unique(genes_by_fold[[i]])
      b <- unique(genes_by_fold[[j]])
      jac[[k]] <- data.frame(
        analysis = label,
        fold_a = fold_names[i],
        fold_b = fold_names[j],
        n_a = length(a),
        n_b = length(b),
        intersection = length(intersect(a, b)),
        union = length(union(a, b)),
        jaccard = ifelse(length(union(a, b)) == 0, NA_real_, length(intersect(a, b)) / length(union(a, b))),
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  list(freq = freq, jaccard = do.call(rbind, jac))
}

bootstrap_ci <- function(scores, n = n_boot) {
  set.seed(20260721)
  young_idx <- which(scores$age_group == "Young")
  old_idx <- which(scores$age_group == "Old")
  boot <- data.frame(
    spearman_all = rep(NA_real_, n),
    spearman_old_only = rep(NA_real_, n),
    young_old_auc = rep(NA_real_, n),
    young_minus_old_median = rep(NA_real_, n)
  )
  for (b in seq_len(n)) {
    idx <- c(sample(young_idx, length(young_idx), replace = TRUE), sample(old_idx, length(old_idx), replace = TRUE))
    s <- scores[idx, ]
    boot$spearman_all[b] <- safe_cor(s$predicted_young_probability, s$age_months)
    old_s <- s[s$age_group == "Old", ]
    boot$spearman_old_only[b] <- safe_cor(old_s$predicted_young_probability, old_s$age_months)
    boot$young_old_auc[b] <- auc_rank(as.integer(s$age_group == "Young"), s$predicted_young_probability)
    boot$young_minus_old_median[b] <- median(s$predicted_young_probability[s$age_group == "Young"]) -
      median(s$predicted_young_probability[s$age_group == "Old"])
  }
  estimates <- c(
    spearman_all = safe_cor(scores$predicted_young_probability, scores$age_months),
    spearman_old_only = safe_cor(
      scores$predicted_young_probability[scores$age_group == "Old"],
      scores$age_months[scores$age_group == "Old"]
    ),
    young_old_auc = auc_rank(as.integer(scores$age_group == "Young"), scores$predicted_young_probability),
    young_minus_old_median =
      median(scores$predicted_young_probability[scores$age_group == "Young"]) -
      median(scores$predicted_young_probability[scores$age_group == "Old"])
  )
  ci <- do.call(rbind, lapply(names(boot), function(metric) {
    vals <- boot[[metric]]
    data.frame(
      metric = metric,
      estimate = unname(estimates[metric]),
      bootstrap_median = median(vals, na.rm = TRUE),
      ci_low = unname(quantile(vals, 0.025, na.rm = TRUE)),
      ci_high = unname(quantile(vals, 0.975, na.rm = TRUE)),
      n_bootstrap = n,
      stringsAsFactors = FALSE
    )
  }))
  list(samples = boot, ci = ci)
}

message("Reading parsed FACS Limb Muscle MSC inputs")
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata$mouse <- metadata$mouse.id
metadata$age_months <- as.integer(sub("m$", "", metadata$age))
metadata$young_label <- as.integer(metadata$age_group == "Young")

message("Opening BPCells matrix")
mat <- open_matrix_dir(matrix_dir)
stopifnot(ncol(mat) == nrow(metadata))
stopifnot(identical(as.character(colnames(mat)), as.character(metadata$index)))

message("Constructing pseudobulk in memory")
mouse_meta <- do.call(rbind, lapply(unique(metadata$mouse), function(mouse_id) {
  idx <- which(metadata$mouse == mouse_id)
  data.frame(
    mouse = mouse_id,
    age = unique(metadata$age[idx]),
    age_months = unique(metadata$age_months[idx]),
    age_group = unique(metadata$age_group[idx]),
    sex = unique(metadata$sex[idx]),
    n_cells = length(idx),
    stringsAsFactors = FALSE
  )
}))
mouse_meta <- mouse_meta[order(mouse_meta$age_months, mouse_meta$sex, mouse_meta$mouse), ]
pseudobulk <- sapply(mouse_meta$mouse, function(mouse_id) {
  idx <- which(metadata$mouse == mouse_id)
  as.numeric(rowSums(mat[, idx, drop = FALSE]))
})
rownames(pseudobulk) <- rownames(mat)
colnames(pseudobulk) <- mouse_meta$mouse
mouse_meta$pseudobulk_library_size <- as.numeric(colSums(pseudobulk))
mouse_meta$pseudobulk_detected_genes <- as.integer(colSums(pseudobulk > 0))
mouse_meta$young_label <- as.integer(mouse_meta$age_group == "Young")

write.csv(mouse_meta, file.path(processed_out, "facs_v2_limb_msc_pseudobulk_metadata.csv"), row.names = FALSE)
saveRDS(pseudobulk, file.path(processed_out, "facs_v2_limb_msc_pseudobulk_counts.rds"))

message("Running primary raw-logCPM nested elastic-net baseline")
primary <- run_nested_en(pseudobulk, mouse_meta, "primary_all_donors", "raw_logcpm")
primary$scores <- primary$scores[order(primary$scores$age_months, primary$scores$sex, primary$scores$mouse), ]

write.csv(primary$scores, file.path(out_root, "facs_v2_elastic_net_heldout_score_points.csv"), row.names = FALSE)
write.csv(primary$tuning, file.path(out_root, "facs_v2_elastic_net_alpha_lambda_by_fold.csv"), row.names = FALSE)
write.csv(primary$coefs, file.path(out_root, "facs_v2_elastic_net_nonzero_coefficients_by_fold.csv"), row.names = FALSE)
write.csv(primary$inner_folds, file.path(out_root, "facs_v2_elastic_net_inner_cv_fold_assignments.csv"), row.names = FALSE)
write.csv(primary$inner_fold_balance, file.path(out_root, "facs_v2_elastic_net_inner_cv_stratification.csv"), row.names = FALSE)

primary_summary <- metric_summary(primary$scores)
write.csv(primary_summary, file.path(out_root, "facs_v2_elastic_net_summary.csv"), row.names = FALSE)

alpha_lambda_summary <- do.call(rbind, lapply(split(primary$scores, primary$scores$analysis), function(s) {
  data.frame(
    analysis = unique(s$analysis),
    selected_alpha_values = paste(sort(unique(s$alpha)), collapse = ";"),
    alpha_min = min(s$alpha),
    alpha_median = median(s$alpha),
    alpha_max = max(s$alpha),
    lambda_used_min = min(s$lambda_used),
    lambda_used_median = median(s$lambda_used),
    lambda_used_max = max(s$lambda_used),
    stringsAsFactors = FALSE
  )
}))
write.csv(alpha_lambda_summary, file.path(out_root, "facs_v2_elastic_net_alpha_lambda_distribution.csv"), row.names = FALSE)

nonzero_summary <- data.frame(
  analysis = "primary_all_donors",
  min_nonzero_genes = min(primary$scores$n_nonzero_genes),
  median_nonzero_genes = median(primary$scores$n_nonzero_genes),
  max_nonzero_genes = max(primary$scores$n_nonzero_genes),
  folds_with_zero_nonzero_genes = sum(primary$scores$n_nonzero_genes == 0),
  stringsAsFactors = FALSE
)
write.csv(nonzero_summary, file.path(out_root, "facs_v2_elastic_net_nonzero_gene_count_summary.csv"), row.names = FALSE)

sel <- selection_stats(primary$coefs, "primary_all_donors")
write.csv(sel$freq, file.path(out_root, "facs_v2_elastic_net_gene_selection_frequency.csv"), row.names = FALSE)
write.csv(sel$jaccard, file.path(out_root, "facs_v2_elastic_net_fold_jaccard.csv"), row.names = FALSE)

library_age <- data.frame(
  metric = c(
    "spearman_library_size_vs_age_all",
    "spearman_library_size_vs_age_old_only",
    "spearman_detected_genes_vs_age_all",
    "spearman_detected_genes_vs_age_old_only",
    "spearman_cell_count_vs_age_all",
    "spearman_cell_count_vs_age_old_only"
  ),
  value = c(
    safe_cor(mouse_meta$pseudobulk_library_size, mouse_meta$age_months),
    safe_cor(mouse_meta$pseudobulk_library_size[mouse_meta$age_group == "Old"], mouse_meta$age_months[mouse_meta$age_group == "Old"]),
    safe_cor(mouse_meta$pseudobulk_detected_genes, mouse_meta$age_months),
    safe_cor(mouse_meta$pseudobulk_detected_genes[mouse_meta$age_group == "Old"], mouse_meta$age_months[mouse_meta$age_group == "Old"]),
    safe_cor(mouse_meta$n_cells, mouse_meta$age_months),
    safe_cor(mouse_meta$n_cells[mouse_meta$age_group == "Old"], mouse_meta$age_months[mouse_meta$age_group == "Old"])
  ),
  stringsAsFactors = FALSE
)
write.csv(library_age, file.path(out_root, "facs_v2_library_size_age_associations.csv"), row.names = FALSE)

library_group <- aggregate(
  cbind(pseudobulk_library_size, pseudobulk_detected_genes, n_cells) ~ age + age_months + age_group + sex,
  data = mouse_meta,
  FUN = function(x) c(n = length(x), median = median(x), min = min(x), max = max(x))
)
library_group <- do.call(data.frame, library_group)
write.csv(library_group, file.path(out_root, "facs_v2_library_size_by_age_group.csv"), row.names = FALSE)

lowest_depth_mouse <- mouse_meta$mouse[which.min(mouse_meta$pseudobulk_library_size)]
message(sprintf("Running low-depth sensitivity without %s", lowest_depth_mouse))
keep_mouse <- mouse_meta$mouse != lowest_depth_mouse
sens <- run_nested_en(pseudobulk[, keep_mouse, drop = FALSE], mouse_meta[keep_mouse, ], "remove_lowest_depth_donor", "raw_logcpm")
sens$scores <- sens$scores[order(sens$scores$age_months, sens$scores$sex, sens$scores$mouse), ]
write.csv(sens$scores, file.path(out_root, "facs_v2_elastic_net_remove_lowest_depth_scores.csv"), row.names = FALSE)
write.csv(metric_summary(sens$scores), file.path(out_root, "facs_v2_elastic_net_remove_lowest_depth_summary.csv"), row.names = FALSE)

message("Running donor bootstrap CIs")
boot <- bootstrap_ci(primary$scores, n_boot)
write.csv(boot$samples, file.path(out_root, "facs_v2_elastic_net_donor_bootstrap_samples.csv"), row.names = FALSE)
write.csv(boot$ci, file.path(out_root, "facs_v2_elastic_net_donor_bootstrap_ci.csv"), row.names = FALSE)

if (has_edgeR) {
  message("Running optional TMM-logCPM sensitivity")
  tmm <- run_nested_en(pseudobulk, mouse_meta, "tmm_logcpm_sensitivity", "tmm_logcpm")
  tmm$scores <- tmm$scores[order(tmm$scores$age_months, tmm$scores$sex, tmm$scores$mouse), ]
  write.csv(tmm$scores, file.path(out_root, "facs_v2_elastic_net_tmm_logcpm_scores.csv"), row.names = FALSE)
  write.csv(metric_summary(tmm$scores), file.path(out_root, "facs_v2_elastic_net_tmm_logcpm_summary.csv"), row.names = FALSE)
  write.csv(tmm$tuning, file.path(out_root, "facs_v2_elastic_net_tmm_logcpm_alpha_lambda_by_fold.csv"), row.names = FALSE)
} else {
  writeLines("edgeR unavailable; TMM-logCPM sensitivity was not run.", file.path(out_root, "facs_v2_elastic_net_tmm_logcpm_not_run.txt"))
}

png(file.path(out_root, "facs_v2_elastic_net_heldout_score_points.png"), width = 1200, height = 850, res = 150)
cols <- ifelse(primary$scores$sex == "female", "#B6465F", "#2878B5")
plot(
  primary$scores$age_months,
  primary$scores$predicted_young_probability,
  pch = ifelse(primary$scores$age_group == "Young", 16, 17),
  col = cols,
  xlab = "Age months",
  ylab = "Nested LOMO elastic-net score",
  main = "FACS v2 Elastic-Net Baseline: Held-out Donor Scores",
  ylim = c(0, 1)
)
text(primary$scores$age_months, primary$scores$predicted_young_probability, labels = primary$scores$mouse, pos = 3, cex = 0.55)
legend("topright", legend = c("female", "male", "Young", "Old"), col = c("#B6465F", "#2878B5", "black", "black"), pch = c(16, 16, 16, 17), bty = "n")
dev.off()

png(file.path(out_root, "facs_v2_library_size_by_age.png"), width = 1200, height = 850, res = 150)
plot(
  mouse_meta$age_months,
  mouse_meta$pseudobulk_library_size,
  pch = ifelse(mouse_meta$sex == "female", 16, 17),
  col = ifelse(mouse_meta$age_group == "Young", "#3A7D44", "#8E3B46"),
  xlab = "Age months",
  ylab = "Pseudobulk library size",
  main = "FACS Limb MSC Pseudobulk Library Size by Age"
)
text(mouse_meta$age_months, mouse_meta$pseudobulk_library_size, labels = mouse_meta$mouse, pos = 3, cex = 0.55)
dev.off()

report_lines <- c(
  "# FACS v2 Elastic-Net Baseline Audit",
  "",
  "This audit reads parsed TMS FACS Limb Muscle MSC data from `data_facs` but writes all generated scripts/outputs outside `data_facs`.",
  "",
  "No GSE176206 data or Kaile frozen model was used.",
  "",
  "## Checks requested",
  "",
  sprintf("1. Inner CV age-stratified: `%s` for all primary folds.", all(primary$scores$inner_cv_age_stratified)),
  "2. Best alpha/lambda by fold: `facs_v2_elastic_net_alpha_lambda_by_fold.csv`; distribution summary: `facs_v2_elastic_net_alpha_lambda_distribution.csv`.",
  sprintf("3. Nonzero genes per fold: min `%s`, median `%s`, max `%s`; zero-gene folds `%s`.", nonzero_summary$min_nonzero_genes, nonzero_summary$median_nonzero_genes, nonzero_summary$max_nonzero_genes, nonzero_summary$folds_with_zero_nonzero_genes),
  "4. Gene selection frequency and fold Jaccard: `facs_v2_elastic_net_gene_selection_frequency.csv`, `facs_v2_elastic_net_fold_jaccard.csv`.",
  "5. Held-out raw 14 points: `facs_v2_elastic_net_heldout_score_points.csv` and `.png`.",
  "6. Library size versus age: `facs_v2_library_size_age_associations.csv`, `facs_v2_library_size_by_age_group.csv`, `facs_v2_library_size_by_age.png`.",
  sprintf("7. Lowest-depth donor removed: `%s`; sensitivity outputs written.", lowest_depth_mouse),
  sprintf("8. Donor bootstrap CIs: `%s` stratified bootstrap replicates in `facs_v2_elastic_net_donor_bootstrap_ci.csv`.", n_boot),
  sprintf("9. TMM-logCPM sensitivity: `%s`.", ifelse(has_edgeR, "completed", "not run because edgeR was unavailable")),
  "",
  "## Primary nested LOMO summary",
  "",
  paste(capture.output(print(primary_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Important limitations",
  "",
  "- Each outer training fold has only 13 mice, so glmnet warns that one binomial class has fewer than 8 observations.",
  "- This is an internal machine-learning baseline, not the final biological Youth Score model.",
  "- The reported `predicted_young_probability` is a classifier output from tiny donor-level folds and should not be treated as a calibrated probability.",
  "- TMM-logCPM sensitivity uses training-only TMM and an external held-out transform against a training reference donor; the held-out donor is not allowed to modify training normalization.",
  "",
  "## Reproducibility",
  "",
  "- Seed: `20260721`",
  "- Session info: `sessionInfo.txt`"
)
writeLines(report_lines, file.path(out_root, "facs_v2_elastic_net_baseline_audit_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))

message("Done")
message(sprintf("Audit report: %s", file.path(out_root, "facs_v2_elastic_net_baseline_audit_report.md")))
