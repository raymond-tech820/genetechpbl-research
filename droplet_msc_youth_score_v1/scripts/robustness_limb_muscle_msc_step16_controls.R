#!/usr/bin/env Rscript

local_lib <- ".Rlibs"
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})

set.seed(20260717)

counts_path <- "data/processed/tms_limb_msc_pseudobulk_counts.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv"
main_signature_path <- "outputs/scores/step12_candidate_signature_versions.csv"
main_scores_path <- "outputs/scores/step12_13_candidate_scores.csv"
reliability_path <- "outputs/stability/gene_reliability_scores.csv"
out_dir <- "outputs/validation"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_permutations <- 20
n_random_gene_sets <- 300
n_weight_shuffles <- 1000
version_primary <- "Medium"
version_comparator <- "Large"
module_n <- c(Medium = 50L, Large = 100L)
w_max <- 3
low_depth_samples <- c("18-F-50", "18-F-51")
sex_linked_genes <- c("Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2")
sex_sensitivity_ratio_threshold <- 2
sex_sensitivity_abs_delta_threshold <- 1
effect_threshold <- 0.5
cpm_threshold <- 1
min_samples_expressed <- 2

spearman_safe <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

weighted_mean_safe <- function(z, w) {
  ok <- is.finite(z) & is.finite(w) & abs(w) > 0
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * z[ok]) / sum(abs(w[ok]))
}

set_factor_levels <- function(metadata) {
  metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
  metadata$age_group <- factor(metadata$age_group, levels = c("Old", "Young"))
  metadata
}

design_full_rank <- function(metadata) {
  metadata <- set_factor_levels(metadata)
  design <- model.matrix(~ sex + age_group, data = metadata)
  qr(design)$rank == ncol(design)
}

filter_genes <- function(count_matrix) {
  dge <- DGEList(counts = count_matrix)
  keep <- rowSums(cpm(dge) > cpm_threshold) >= min_samples_expressed
  keep[is.na(keep)] <- FALSE
  keep
}

fit_age_model <- function(count_matrix, metadata) {
  metadata <- set_factor_levels(metadata)
  design <- model.matrix(~ sex + age_group, data = metadata)
  if (qr(design)$rank != ncol(design)) {
    return(list(ok = FALSE, reason = "design_not_full_rank"))
  }
  dge <- DGEList(counts = count_matrix, samples = metadata)
  dge <- normLibSizes(dge, method = "TMM")
  fit_result <- tryCatch({
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design, robust = TRUE)
    test <- glmQLFTest(fit, coef = "age_groupYoung")
    table <- topTags(test, n = Inf, sort.by = "none")$table
    table$gene <- rownames(table)
    logcpm <- cpm(dge, log = TRUE, prior.count = 2)
    age_rho <- apply(logcpm, 1, function(x) {
      suppressWarnings(cor(x, metadata$age_months, method = "spearman"))
    })
    list(ok = TRUE, dge = dge, table = table, logcpm = logcpm, age_rho = age_rho, design = design)
  }, error = function(e) {
    list(ok = FALSE, reason = conditionMessage(e))
  })
  fit_result
}

build_reliability <- function(count_matrix, metadata, no_sex_filter = FALSE) {
  keep <- filter_genes(count_matrix)
  counts_f <- count_matrix[keep, metadata$sample_id, drop = FALSE]
  fit <- fit_age_model(counts_f, metadata)
  if (!isTRUE(fit$ok)) {
    return(list(ok = FALSE, reason = fit$reason))
  }

  genes <- rownames(counts_f)
  de <- fit$table[match(genes, fit$table$gene), ]
  logfc <- de$logFC
  names(logfc) <- genes
  fdr <- de$FDR
  names(fdr) <- genes
  age_rho <- fit$age_rho[genes]

  lomo_logfc <- matrix(NA_real_, nrow = length(genes), ncol = nrow(metadata),
                       dimnames = list(genes, metadata$sample_id))
  lomo_valid <- setNames(rep(FALSE, nrow(metadata)), metadata$sample_id)
  for (sample_id in metadata$sample_id) {
    inner_meta <- metadata[metadata$sample_id != sample_id, , drop = FALSE]
    if (!design_full_rank(inner_meta)) next
    inner_fit <- fit_age_model(counts_f[, inner_meta$sample_id, drop = FALSE], inner_meta)
    if (!isTRUE(inner_fit$ok)) next
    lomo_valid[sample_id] <- TRUE
    lomo_logfc[inner_fit$table$gene, sample_id] <- inner_fit$table$logFC
  }
  full_sign <- sign(logfc[genes])
  sign_match <- sweep(sign(lomo_logfc), 1, full_sign, FUN = "==")
  sign_match[full_sign == 0, ] <- TRUE
  sign_match[, !lomo_valid] <- NA
  lomo_sign_rate <- rowMeans(sign_match, na.rm = TRUE)
  valid_lomo_count <- rowSums(!is.na(sign_match))
  lomo_sign_rate[valid_lomo_count == 0] <- NA_real_

  young_ids <- metadata$sample_id[metadata$age_group == "Young"]
  valid_young_ids <- young_ids[lomo_valid[young_ids]]
  if (length(valid_young_ids) > 0) {
    young_reversal <- apply(!sign_match[, valid_young_ids, drop = FALSE], 1, any, na.rm = TRUE)
  } else {
    young_reversal <- rep(FALSE, length(genes))
  }

  low_depth_removed <- intersect(low_depth_samples, metadata$sample_id)
  low_depth_triggered <- length(low_depth_removed) > 0
  low_depth_valid <- FALSE
  low_depth_sign_match <- setNames(rep(TRUE, length(genes)), genes)
  low_depth_abs_change <- setNames(rep(NA_real_, length(genes)), genes)
  if (low_depth_triggered) {
    reduced_meta <- metadata[!(metadata$sample_id %in% low_depth_removed), , drop = FALSE]
    if (design_full_rank(reduced_meta)) {
      reduced_fit <- fit_age_model(counts_f[, reduced_meta$sample_id, drop = FALSE], reduced_meta)
      if (isTRUE(reduced_fit$ok)) {
        low_depth_valid <- TRUE
        reduced_logfc <- reduced_fit$table$logFC
        names(reduced_logfc) <- reduced_fit$table$gene
        low_depth_sign_match[genes] <- sign(logfc[genes]) == sign(reduced_logfc[genes])
        low_depth_abs_change[genes] <- abs(reduced_logfc[genes] - logfc[genes])
      }
    }
  }

  old_male <- metadata$age_group == "Old" & metadata$sex == "male"
  old_female <- metadata$age_group == "Old" & metadata$sex == "female"
  if (sum(old_male) > 0 && sum(old_female) > 0) {
    old_sex_delta <- rowMeans(fit$logcpm[, old_male, drop = FALSE]) -
      rowMeans(fit$logcpm[, old_female, drop = FALSE])
  } else {
    old_sex_delta <- setNames(rep(NA_real_, length(genes)), genes)
  }
  sex_ratio <- abs(old_sex_delta[genes]) / pmax(abs(logfc[genes]), 1e-6)
  strong_sex_sensitive <- abs(old_sex_delta[genes]) > sex_sensitivity_abs_delta_threshold &
    sex_ratio > sex_sensitivity_ratio_threshold
  strong_sex_sensitive[is.na(strong_sex_sensitive)] <- FALSE
  sex_linked <- genes %in% sex_linked_genes

  trend_compatible <- ifelse(logfc[genes] > 0, age_rho[genes] < 0, age_rho[genes] > 0)
  trend_compatible[is.na(trend_compatible)] <- FALSE

  sex_pass <- if (no_sex_filter) rep(TRUE, length(genes)) else (!sex_linked & !strong_sex_sensitive)
  pi_sex <- if (no_sex_filter) rep(1, length(genes)) else ifelse(sex_pass, 1, 0)
  lomo_pass <- !is.na(lomo_sign_rate) & lomo_sign_rate >= 0.9 & !young_reversal
  low_depth_pass <- as.logical(low_depth_sign_match[genes])
  low_depth_pass[is.na(low_depth_pass)] <- FALSE

  reliability <- data.frame(
    gene = genes,
    adjusted_logFC = logfc[genes],
    FDR = fdr[genes],
    age_rho = age_rho[genes],
    direction = ifelse(logfc[genes] > 0, "young_high", "old_high"),
    continuous_age_trend_compatible = trend_compatible,
    LOMO_sign_rate = lomo_sign_rate,
    valid_inner_LOMO_count = valid_lomo_count,
    young_mouse_removal_direction_reversal = young_reversal,
    low_depth_sign_match = low_depth_pass,
    low_depth_abs_effect_size_change = low_depth_abs_change[genes],
    low_depth_triggered = low_depth_triggered,
    low_depth_valid = low_depth_valid,
    sex_effect_old_male_minus_old_female_logcpm = old_sex_delta[genes],
    sex_effect_ratio_vs_age = sex_ratio,
    obvious_sex_linked_gene = sex_linked,
    strong_sex_sensitive_gene = strong_sex_sensitive,
    effect_size_pass_abs_logFC_gt_0_5 = abs(logfc[genes]) > effect_threshold,
    LOMO_pass = lomo_pass,
    low_depth_pass = low_depth_pass,
    sex_pass = sex_pass,
    pi_LOMO = pmax(0, pmin(1, lomo_sign_rate)),
    pi_depth = ifelse(low_depth_pass, 1, 0),
    pi_sex = pi_sex,
    stringsAsFactors = FALSE
  )
  reliability$r_g <- reliability$pi_LOMO * reliability$pi_depth * reliability$pi_sex
  reliability$passes_step9_initial_reliability <- reliability$effect_size_pass_abs_logFC_gt_0_5 &
    reliability$continuous_age_trend_compatible &
    reliability$LOMO_pass &
    reliability$low_depth_pass &
    reliability$sex_pass
  reliability$q_g <- abs(reliability$adjusted_logFC) * abs(reliability$age_rho) * reliability$r_g

  list(ok = TRUE, counts = counts_f, fit = fit, reliability = reliability)
}

select_signature <- function(reliability, version, n_per_direction) {
  candidate <- reliability[reliability$passes_step9_initial_reliability, ]
  young <- candidate[candidate$direction == "young_high", ]
  old <- candidate[candidate$direction == "old_high", ]
  young <- young[order(-young$q_g, -abs(young$adjusted_logFC), -abs(young$age_rho), young$gene), ]
  old <- old[order(-old$q_g, -abs(old$adjusted_logFC), -abs(old$age_rho), old$gene), ]
  sig <- rbind(head(young, n_per_direction), head(old, n_per_direction))
  sig$version <- version
  sig$module <- ifelse(sig$direction == "young_high", "young_module", "old_module")
  sig$weight <- pmin(abs(sig$adjusted_logFC), w_max) * sig$r_g
  sig
}

score_signature <- function(signature, logcpm, metadata, calibration_metadata = metadata) {
  sig <- signature[signature$gene %in% rownames(logcpm), ]
  sig$training_mean <- rowMeans(logcpm[sig$gene, , drop = FALSE])
  sig$training_sd <- apply(logcpm[sig$gene, , drop = FALSE], 1, sd)
  sig <- sig[is.finite(sig$training_sd) & sig$training_sd > 0 & is.finite(sig$weight) & sig$weight > 0, ]
  young_sig <- sig[sig$module == "young_module", ]
  old_sig <- sig[sig$module == "old_module", ]
  raw <- setNames(rep(NA_real_, ncol(logcpm)), colnames(logcpm))
  for (sample_id in colnames(logcpm)) {
    yz <- (as.numeric(logcpm[young_sig$gene, sample_id]) - young_sig$training_mean) / young_sig$training_sd
    oz <- (as.numeric(logcpm[old_sig$gene, sample_id]) - old_sig$training_mean) / old_sig$training_sd
    raw[sample_id] <- weighted_mean_safe(yz, young_sig$weight) - weighted_mean_safe(oz, old_sig$weight)
  }
  age_group <- calibration_metadata$age_group[match(names(raw), calibration_metadata$sample_id)]
  young_center <- median(raw[age_group == "Young"], na.rm = TRUE)
  old_center <- median(raw[age_group == "Old"], na.rm = TRUE)
  denom <- young_center - old_center
  calibrated <- (raw - old_center) / denom
  data.frame(
    sample_id = names(raw),
    raw_score = raw,
    calibrated_score = calibrated,
    young_reference_center = young_center,
    old_reference_center = old_center,
    calibration_denominator = denom,
    stringsAsFactors = FALSE
  )
}

score_metric <- function(score_df, truth_metadata) {
  x <- merge(score_df, truth_metadata, by.x = "sample_id", by.y = "sample_id", all.x = TRUE)
  young <- x$calibrated_score[x$age_group == "Young"]
  old <- x$calibrated_score[x$age_group == "Old"]
  raw_delta <- if ("raw_score" %in% colnames(x)) {
    median(x$raw_score[x$age_group == "Young"], na.rm = TRUE) -
      median(x$raw_score[x$age_group == "Old"], na.rm = TRUE)
  } else {
    NA_real_
  }
  data.frame(
    rho_age = spearman_safe(x$calibrated_score, x$age_months),
    abs_rho_age = abs(spearman_safe(x$calibrated_score, x$age_months)),
    young_old_delta = median(young, na.rm = TRUE) - median(old, na.rm = TRUE),
    raw_young_old_delta = raw_delta,
    old_male_median = median(x$calibrated_score[x$age_group == "Old" & x$sex == "male"], na.rm = TRUE),
    old_female_median = median(x$calibrated_score[x$age_group == "Old" & x$sex == "female"], na.rm = TRUE),
    old_male_minus_female = median(x$calibrated_score[x$age_group == "Old" & x$sex == "male"], na.rm = TRUE) -
      median(x$calibrated_score[x$age_group == "Old" & x$sex == "female"], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

bind_rows_fill <- function(frames) {
  all_cols <- unique(unlist(lapply(frames, names)))
  padded <- lapply(frames, function(x) {
    missing <- setdiff(all_cols, names(x))
    for (col in missing) x[[col]] <- NA
    x[, all_cols, drop = FALSE]
  })
  do.call(rbind, padded)
}

score_all_samples_with_train_fit <- function(signature, train_result, train_metadata, all_counts, all_metadata) {
  genes <- rownames(train_result$counts)
  train_ids <- train_metadata$sample_id
  train_lib <- train_result$fit$dge$samples$lib.size * train_result$fit$dge$samples$norm.factors
  names(train_lib) <- train_ids
  ref_sample <- names(train_lib)[which.min(abs(train_lib - median(train_lib)))]
  combined_ids <- all_metadata$sample_id
  combined_counts <- all_counts[genes, combined_ids, drop = FALSE]
  dge <- DGEList(counts = combined_counts)
  dge <- normLibSizes(dge, method = "TMM", refColumn = match(ref_sample, combined_ids))
  logcpm <- cpm(dge, log = TRUE, prior.count = 2)
  score_signature(signature, logcpm, all_metadata, calibration_metadata = train_metadata)
}

message("Reading inputs")
counts <- readRDS(counts_path)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
main_signature <- read.csv(main_signature_path, check.names = FALSE, stringsAsFactors = FALSE)
main_scores <- read.csv(main_scores_path, check.names = FALSE, stringsAsFactors = FALSE)
full_reliability <- read.csv(reliability_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!identical(colnames(counts), metadata$sample_id)) {
  stop("Counts columns do not match metadata sample_id order")
}

main_medium_scores <- main_scores[main_scores$version == version_primary, ]
main_medium_score_df <- data.frame(
  sample_id = main_medium_scores$sample_id,
  calibrated_score = main_medium_scores$youth_score_raw_calibrated,
  raw_score = main_medium_scores$score_raw,
  stringsAsFactors = FALSE
)
main_medium_metric <- score_metric(main_medium_score_df, metadata)

message("Step 16.1: age-label permutation controls")
permutation_file <- file.path(out_dir, "step16_age_label_permutation_null.csv")
if (file.exists(permutation_file)) {
  permutation <- read.csv(permutation_file, check.names = FALSE, stringsAsFactors = FALSE)
  if ("raw_young_old_delta" %in% colnames(permutation)) {
    message("  reusing existing permutation output")
  } else {
    message("  existing permutation output lacks raw delta; recomputing")
    file.remove(permutation_file)
    permutation <- NULL
  }
}
if (!exists("permutation") || is.null(permutation)) {
  permutation_rows <- list()
  for (b in seq_len(n_permutations)) {
    message(sprintf("  permutation %s / %s", b, n_permutations))
    perm_metadata <- metadata
    perm_idx <- sample(seq_len(nrow(metadata)))
    perm_metadata$age_months <- metadata$age_months[perm_idx]
    perm_metadata$age <- metadata$age[perm_idx]
    perm_metadata$age_group <- metadata$age_group[perm_idx]
    perm_metadata$young_label <- ifelse(perm_metadata$age_group == "Young", 1, 0)
    perm_metadata$old_label <- ifelse(perm_metadata$age_group == "Old", 1, 0)
    perm_metadata$age_group_model <- perm_metadata$age_group
    perm_metadata$age_months_continuous <- perm_metadata$age_months

    if (!design_full_rank(perm_metadata)) {
      permutation_rows[[length(permutation_rows) + 1]] <- data.frame(
        permutation = b,
        version = NA_character_,
        ok = FALSE,
        reason = "permuted_design_not_full_rank",
        stringsAsFactors = FALSE
      )
      next
    }
    perm_result <- build_reliability(counts, perm_metadata, no_sex_filter = FALSE)
    if (!isTRUE(perm_result$ok)) {
      permutation_rows[[length(permutation_rows) + 1]] <- data.frame(
        permutation = b,
        version = NA_character_,
        ok = FALSE,
        reason = perm_result$reason,
        stringsAsFactors = FALSE
      )
      next
    }
    for (version_name in c(version_primary, version_comparator)) {
      sig <- select_signature(perm_result$reliability, version_name, module_n[[version_name]])
      if (sum(sig$direction == "young_high") < module_n[[version_name]] ||
          sum(sig$direction == "old_high") < module_n[[version_name]]) {
        permutation_rows[[length(permutation_rows) + 1]] <- data.frame(
          permutation = b,
          version = version_name,
          ok = FALSE,
          reason = "not_enough_candidates",
          stringsAsFactors = FALSE
        )
        next
      }
      score <- score_signature(sig, perm_result$fit$logcpm, perm_metadata, calibration_metadata = perm_metadata)
      metric <- score_metric(score, metadata)
      permutation_rows[[length(permutation_rows) + 1]] <- data.frame(
        permutation = b,
        version = version_name,
        ok = TRUE,
        reason = "ok",
        n_candidate_genes = sum(perm_result$reliability$passes_step9_initial_reliability),
        selected_young = sum(sig$direction == "young_high"),
        selected_old = sum(sig$direction == "old_high"),
        metric,
        stringsAsFactors = FALSE
      )
    }
  }
  permutation <- bind_rows_fill(permutation_rows)
  write.csv(permutation, permutation_file, row.names = FALSE)
}

message("Step 16.2: expression-matched random gene-set controls")
medium_sig <- main_signature[main_signature$version == version_primary, ]
medium_sig$direction <- ifelse(medium_sig$module == "young_module", "young_high", "old_high")
pipeline_main <- build_reliability(counts, metadata, no_sex_filter = FALSE)
if (!isTRUE(pipeline_main$ok)) stop(pipeline_main$reason)
logcpm <- pipeline_main$fit$logcpm
mean_expr <- rowMeans(logcpm)
expr_bins <- cut(mean_expr, breaks = quantile(mean_expr, probs = seq(0, 1, length.out = 11), na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE)
names(expr_bins) <- names(mean_expr)

sample_matched_genes <- function(target_genes, exclude_genes) {
  picked <- character(0)
  for (gene in target_genes) {
    bin <- expr_bins[gene]
    pool <- names(expr_bins)[expr_bins == bin & !(names(expr_bins) %in% c(exclude_genes, picked))]
    if (length(pool) == 0) {
      pool <- names(expr_bins)[!(names(expr_bins) %in% c(exclude_genes, picked))]
    }
    picked <- c(picked, sample(pool, 1))
  }
  picked
}

random_file <- file.path(out_dir, "step16_random_gene_set_controls.csv")
if (file.exists(random_file)) {
  random_controls <- read.csv(random_file, check.names = FALSE, stringsAsFactors = FALSE)
  if ("raw_young_old_delta" %in% colnames(random_controls)) {
    message("  reusing existing random gene-set output")
  } else {
    message("  existing random gene-set output lacks raw delta; recomputing")
    file.remove(random_file)
    random_controls <- NULL
  }
}
if (!exists("random_controls") || is.null(random_controls)) {
  random_rows <- list()
  medium_weights_y <- medium_sig$weight[medium_sig$module == "young_module"]
  medium_weights_o <- medium_sig$weight[medium_sig$module == "old_module"]
  for (b in seq_len(n_random_gene_sets)) {
    rand_y <- sample_matched_genes(medium_sig$gene[medium_sig$module == "young_module"], medium_sig$gene)
    rand_o <- sample_matched_genes(medium_sig$gene[medium_sig$module == "old_module"], c(medium_sig$gene, rand_y))
    sig <- data.frame(
      gene = c(rand_y, rand_o),
      module = rep(c("young_module", "old_module"), each = module_n[[version_primary]]),
      direction = rep(c("young_high", "old_high"), each = module_n[[version_primary]]),
      weight = c(sample(medium_weights_y), sample(medium_weights_o)),
      stringsAsFactors = FALSE
    )
    score <- score_signature(sig, logcpm, metadata, calibration_metadata = metadata)
    metric <- score_metric(score, metadata)
    random_rows[[b]] <- data.frame(control_id = b, metric, stringsAsFactors = FALSE)
  }
  random_controls <- do.call(rbind, random_rows)
  write.csv(random_controls, random_file, row.names = FALSE)
}

message("Step 16.3: Medium weight-shuffle controls")
weight_file <- file.path(out_dir, "step16_weight_shuffle_controls.csv")
if (file.exists(weight_file)) {
  weight_shuffle <- read.csv(weight_file, check.names = FALSE, stringsAsFactors = FALSE)
  if ("raw_young_old_delta" %in% colnames(weight_shuffle)) {
    message("  reusing existing weight-shuffle output")
  } else {
    message("  existing weight-shuffle output lacks raw delta; recomputing")
    file.remove(weight_file)
    weight_shuffle <- NULL
  }
}
if (!exists("weight_shuffle") || is.null(weight_shuffle)) {
  weight_rows <- list()
  for (b in seq_len(n_weight_shuffles)) {
    sig <- medium_sig
    sig$weight[sig$module == "young_module"] <- sample(sig$weight[sig$module == "young_module"])
    sig$weight[sig$module == "old_module"] <- sample(sig$weight[sig$module == "old_module"])
    score <- score_signature(sig, logcpm, metadata, calibration_metadata = metadata)
    metric <- score_metric(score, metadata)
    weight_rows[[b]] <- data.frame(control_id = b, metric, stringsAsFactors = FALSE)
  }
  weight_shuffle <- do.call(rbind, weight_rows)
  write.csv(weight_shuffle, weight_file, row.names = FALSE)
}

message("Step 16.4: low-depth score sensitivity")
reduced_metadata <- metadata[!(metadata$sample_id %in% low_depth_samples), , drop = FALSE]
low_depth_result <- build_reliability(counts[, reduced_metadata$sample_id, drop = FALSE], reduced_metadata, no_sex_filter = FALSE)
if (!isTRUE(low_depth_result$ok)) stop(low_depth_result$reason)
low_depth_rows <- list()
for (version_name in c(version_primary, version_comparator)) {
  sig <- select_signature(low_depth_result$reliability, version_name, module_n[[version_name]])
  score <- score_all_samples_with_train_fit(sig, low_depth_result, reduced_metadata, counts, metadata)
  metric <- score_metric(score, metadata)
  main_version_scores <- main_scores[main_scores$version == version_name, c("sample_id", "youth_score_raw_calibrated")]
  names(main_version_scores)[2] <- "main_full_score"
  merged <- merge(main_version_scores, score[, c("sample_id", "calibrated_score")], by = "sample_id")
  common <- merged[!(merged$sample_id %in% low_depth_samples), ]
  low_depth_rows[[version_name]] <- data.frame(
    version = version_name,
    reduced_signature_size = nrow(sig),
    reduced_selected_young = sum(sig$direction == "young_high"),
    reduced_selected_old = sum(sig$direction == "old_high"),
    score_spearman_all_mice = spearman_safe(merged$main_full_score, merged$calibrated_score),
    score_spearman_excluding_low_depth_mice = spearman_safe(common$main_full_score, common$calibrated_score),
    rank_spearman_all_mice = spearman_safe(rank(merged$main_full_score), rank(merged$calibrated_score)),
    rank_spearman_excluding_low_depth_mice = spearman_safe(rank(common$main_full_score), rank(common$calibrated_score)),
    metric,
    stringsAsFactors = FALSE
  )
}
low_depth_sensitivity <- do.call(rbind, low_depth_rows)
write.csv(low_depth_sensitivity, file.path(out_dir, "step16_low_depth_score_sensitivity.csv"), row.names = FALSE)

message("Step 16.5: sex-linked inclusion stress test")
sex_stress_result <- build_reliability(counts, metadata, no_sex_filter = TRUE)
if (!isTRUE(sex_stress_result$ok)) stop(sex_stress_result$reason)
sex_rows <- list()
for (version_name in c(version_primary, version_comparator)) {
  sig <- select_signature(sex_stress_result$reliability, version_name, module_n[[version_name]])
  score <- score_signature(sig, sex_stress_result$fit$logcpm, metadata, calibration_metadata = metadata)
  metric <- score_metric(score, metadata)
  sex_rows[[version_name]] <- data.frame(
    version = version_name,
    model = "no_sex_filter_stress",
    selected_signature_size = nrow(sig),
    selected_obvious_sex_linked = sum(sig$obvious_sex_linked_gene),
    selected_strong_sex_sensitive = sum(sig$strong_sex_sensitive_gene),
    metric,
    stringsAsFactors = FALSE
  )
}
main_sex_rows <- lapply(c(version_primary, version_comparator), function(version_name) {
  x <- main_scores[main_scores$version == version_name, ]
  score <- data.frame(
    sample_id = x$sample_id,
    calibrated_score = x$youth_score_raw_calibrated,
    raw_score = x$score_raw,
    stringsAsFactors = FALSE
  )
  metric <- score_metric(score, metadata)
  data.frame(
    version = version_name,
    model = "main_sex_filtered",
    selected_signature_size = nrow(main_signature[main_signature$version == version_name, ]),
    selected_obvious_sex_linked = sum(main_signature$version == version_name & main_signature$gene %in% sex_linked_genes),
    selected_strong_sex_sensitive = NA_integer_,
    metric,
    stringsAsFactors = FALSE
  )
})
sex_robustness <- do.call(rbind, c(main_sex_rows, sex_rows))
write.csv(sex_robustness, file.path(out_dir, "step16_sex_robustness_stress.csv"), row.names = FALSE)

message("Summarizing controls")
permutation_summary <- do.call(rbind, lapply(c(version_primary, version_comparator), function(version_name) {
  x <- permutation[permutation$version == version_name & permutation$ok, ]
  true_score <- main_scores[main_scores$version == version_name, ]
  true_metric <- score_metric(data.frame(
    sample_id = true_score$sample_id,
    calibrated_score = true_score$youth_score_raw_calibrated,
    raw_score = true_score$score_raw,
    stringsAsFactors = FALSE
  ), metadata)
  data.frame(
    version = version_name,
    control = "age_label_permutation_training_pipeline",
    n_successful_permutations = nrow(x),
    true_abs_rho_age = true_metric$abs_rho_age,
    null_abs_rho_age_median = median(x$abs_rho_age, na.rm = TRUE),
    null_abs_rho_age_95pct = quantile(x$abs_rho_age, 0.95, na.rm = TRUE),
    empirical_p_abs_rho = (sum(x$abs_rho_age >= true_metric$abs_rho_age, na.rm = TRUE) + 1) / (nrow(x) + 1),
    true_young_old_delta_calibrated = true_metric$young_old_delta,
    true_raw_young_old_delta = true_metric$raw_young_old_delta,
    null_raw_delta_median = median(x$raw_young_old_delta, na.rm = TRUE),
    null_raw_delta_95pct = quantile(x$raw_young_old_delta, 0.95, na.rm = TRUE),
    empirical_p_raw_delta = (sum(x$raw_young_old_delta >= true_metric$raw_young_old_delta, na.rm = TRUE) + 1) / (nrow(x) + 1),
    note = "reran training-side feature selection; not nested outer-LOMO null",
    stringsAsFactors = FALSE
  )
}))

random_summary <- data.frame(
  version = version_primary,
  control = "expression_matched_random_gene_sets",
  n_controls = nrow(random_controls),
  true_abs_rho_age = main_medium_metric$abs_rho_age,
  null_abs_rho_age_median = median(random_controls$abs_rho_age, na.rm = TRUE),
  null_abs_rho_age_95pct = quantile(random_controls$abs_rho_age, 0.95, na.rm = TRUE),
  empirical_p_abs_rho = (sum(random_controls$abs_rho_age >= main_medium_metric$abs_rho_age, na.rm = TRUE) + 1) /
    (nrow(random_controls) + 1),
  true_young_old_delta_calibrated = main_medium_metric$young_old_delta,
  true_raw_young_old_delta = main_medium_metric$raw_young_old_delta,
  null_raw_delta_median = median(random_controls$raw_young_old_delta, na.rm = TRUE),
  null_raw_delta_95pct = quantile(random_controls$raw_young_old_delta, 0.95, na.rm = TRUE),
  empirical_p_raw_delta = (sum(random_controls$raw_young_old_delta >= main_medium_metric$raw_young_old_delta, na.rm = TRUE) + 1) /
    (nrow(random_controls) + 1),
  note = "score-level control matched by expression bins and signature size",
  stringsAsFactors = FALSE
)

weight_summary <- data.frame(
  version = version_primary,
  control = "medium_weight_shuffle",
  n_controls = nrow(weight_shuffle),
  true_abs_rho_age = main_medium_metric$abs_rho_age,
  null_abs_rho_age_median = median(weight_shuffle$abs_rho_age, na.rm = TRUE),
  null_abs_rho_age_95pct = quantile(weight_shuffle$abs_rho_age, 0.95, na.rm = TRUE),
  empirical_p_abs_rho = (sum(weight_shuffle$abs_rho_age >= main_medium_metric$abs_rho_age, na.rm = TRUE) + 1) /
    (nrow(weight_shuffle) + 1),
  true_young_old_delta_calibrated = main_medium_metric$young_old_delta,
  true_raw_young_old_delta = main_medium_metric$raw_young_old_delta,
  null_raw_delta_median = median(weight_shuffle$raw_young_old_delta, na.rm = TRUE),
  null_raw_delta_95pct = quantile(weight_shuffle$raw_young_old_delta, 0.95, na.rm = TRUE),
  empirical_p_raw_delta = (sum(weight_shuffle$raw_young_old_delta >= main_medium_metric$raw_young_old_delta, na.rm = TRUE) + 1) /
    (nrow(weight_shuffle) + 1),
  note = "same Medium genes; weights shuffled within modules",
  stringsAsFactors = FALSE
)

control_summary <- bind_rows_fill(list(permutation_summary, random_summary, weight_summary))
write.csv(control_summary, file.path(out_dir, "step16_control_summary.csv"), row.names = FALSE)

p_null <- ggplot(permutation[permutation$ok & permutation$version %in% c(version_primary, version_comparator), ],
                 aes(x = abs_rho_age, fill = version)) +
  geom_histogram(alpha = 0.65, bins = 20, position = "identity") +
  geom_vline(data = control_summary[control_summary$control == "age_label_permutation_training_pipeline", ],
             aes(xintercept = true_abs_rho_age, color = version), linewidth = 0.8) +
  labs(
    title = "Step 16 age-label permutation null",
    x = "|Spearman rho(score, true age)|",
    y = "Permutation count",
    fill = "Version",
    color = "True"
  ) +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "step16_age_label_permutation_null.png"), p_null, width = 8, height = 5, dpi = 180)

p_controls <- rbind(
  data.frame(control = "random_gene_set", abs_rho_age = random_controls$abs_rho_age),
  data.frame(control = "weight_shuffle", abs_rho_age = weight_shuffle$abs_rho_age)
)
p_ctrl <- ggplot(p_controls, aes(x = abs_rho_age, fill = control)) +
  geom_histogram(alpha = 0.65, bins = 30, position = "identity") +
  geom_vline(xintercept = main_medium_metric$abs_rho_age, color = "black", linewidth = 0.8) +
  labs(
    title = "Step 16 Medium score-level controls",
    x = "|Spearman rho(score, age)|",
    y = "Control count",
    fill = "Control"
  ) +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "step16_score_level_controls.png"), p_ctrl, width = 8, height = 5, dpi = 180)

report_lines <- c(
  "# Step 16: Null and Robustness Controls",
  "",
  "## Scope",
  "",
  "Medium is treated as the primary model and Large as the high-stability comparator. Small is not carried forward except through prior Step 15 documentation.",
  "",
  "## Important Limitation",
  "",
  sprintf("Age-label permutation used %s successful-attempt target iterations and reran training-side TMM, DE, continuous age trend, LOMO stability, low-depth sensitivity, sex-sensitive filtering, ranking, signature selection, and score construction. It is not a full nested outer-LOMO null, so it should be treated as a practical null control rather than the final exhaustive permutation analysis.", n_permutations),
  "",
  "## Control Summary",
  "",
  paste(capture.output(print(control_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Low-Depth Score Sensitivity",
  "",
  paste(capture.output(print(low_depth_sensitivity, row.names = FALSE)), collapse = "\n"),
  "",
  "## Sex Robustness",
  "",
  paste(capture.output(print(sex_robustness, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "- Age-label permutation tests whether the pipeline can recover similar age association after mouse-level age labels are broken.",
  "- Random gene-set controls test whether expression-matched arbitrary genes of the same size can perform similarly.",
  "- Weight shuffling tests whether Medium performance depends on precise weights or mostly on gene identity/module direction.",
  "- Calibrated young-old delta is fixed by training-set median calibration in full-data score-level controls, so the summary uses raw-score young-old delta for the delta null comparisons.",
  "- Low-depth sensitivity compares final score/rank behavior after removing `18-F-50` and `18-F-51` from training.",
  "- Sex robustness checks whether allowing sex-linked/sex-sensitive genes worsens old male-female separation.",
  "",
  "## Outputs",
  "",
  "- `outputs/validation/step16_control_summary.csv`",
  "- `outputs/validation/step16_age_label_permutation_null.csv`",
  "- `outputs/validation/step16_random_gene_set_controls.csv`",
  "- `outputs/validation/step16_weight_shuffle_controls.csv`",
  "- `outputs/validation/step16_low_depth_score_sensitivity.csv`",
  "- `outputs/validation/step16_sex_robustness_stress.csv`",
  "- `outputs/validation/step16_age_label_permutation_null.png`",
  "- `outputs/validation/step16_score_level_controls.png`"
)
writeLines(report_lines, file.path(out_dir, "step16_null_and_robustness_report.md"))

cat("Step 16 controls complete\n")
print(control_summary, row.names = FALSE)
cat("\nLow-depth sensitivity:\n")
print(low_depth_sensitivity, row.names = FALSE)
cat("\nSex robustness:\n")
print(sex_robustness, row.names = FALSE)
