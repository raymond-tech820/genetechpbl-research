#!/usr/bin/env Rscript

local_lib <- ".Rlibs"
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(edgeR)
})

counts_path <- "data/processed/tms_limb_msc_pseudobulk_counts.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv"
validation_dir <- "outputs/validation"

dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

versions <- data.frame(
  version = c("Small", "Medium", "Large"),
  genes_per_direction = c(20L, 50L, 100L),
  stringsAsFactors = FALSE
)

cpm_threshold <- 1
min_samples_expressed <- 2
effect_threshold <- 0.5
w_max <- 3
low_depth_samples <- c("18-F-50", "18-F-51")
sex_linked_genes <- c("Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2")
sex_sensitivity_ratio_threshold <- 2
sex_sensitivity_abs_delta_threshold <- 1

set_factor_levels <- function(metadata) {
  metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
  metadata$age_group <- factor(metadata$age_group, levels = c("Old", "Young"))
  metadata
}

design_for <- function(metadata) {
  metadata <- set_factor_levels(metadata)
  model.matrix(~ sex + age_group, data = metadata)
}

is_design_full_rank <- function(metadata) {
  design <- design_for(metadata)
  qr(design)$rank == ncol(design)
}

filter_genes_train <- function(count_matrix) {
  dge <- DGEList(counts = count_matrix)
  keep <- rowSums(cpm(dge) > cpm_threshold) >= min_samples_expressed
  keep[is.na(keep)] <- FALSE
  keep
}

fit_age_model <- function(count_matrix, sample_metadata) {
  sample_metadata <- set_factor_levels(sample_metadata)
  design <- design_for(sample_metadata)
  if (qr(design)$rank != ncol(design)) {
    return(list(ok = FALSE, reason = "design_not_full_rank", design = design))
  }
  if (!all(c("Old", "Young") %in% as.character(sample_metadata$age_group))) {
    return(list(ok = FALSE, reason = "missing_age_group", design = design))
  }

  dge <- DGEList(counts = count_matrix, samples = sample_metadata)
  dge <- calcNormFactors(dge, method = "TMM")
  fit_result <- tryCatch({
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design, robust = TRUE)
    test <- glmQLFTest(fit, coef = "age_groupYoung")
    table <- topTags(test, n = Inf, sort.by = "none")$table
    table$gene <- rownames(table)
    logcpm <- cpm(dge, log = TRUE, prior.count = 2)
    age_rho <- apply(logcpm, 1, function(x) {
      suppressWarnings(cor(x, sample_metadata$age_months, method = "spearman"))
    })
    list(ok = TRUE, table = table, logcpm = logcpm, age_rho = age_rho, design = design, dge = dge)
  }, error = function(e) {
    list(ok = FALSE, reason = conditionMessage(e), design = design)
  })
  fit_result
}

weighted_mean_safe <- function(z, w) {
  ok <- is.finite(z) & is.finite(w) & abs(w) > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(w[ok] * z[ok]) / sum(abs(w[ok]))
}

score_one_sample <- function(signature, train_logcpm, heldout_logcpm, heldout_sample_id) {
  signature$training_mean <- NA_real_
  signature$training_sd <- NA_real_
  signature$available_in_heldout <- signature$gene %in% rownames(heldout_logcpm)

  for (i in seq_len(nrow(signature))) {
    gene <- signature$gene[i]
    if (gene %in% rownames(train_logcpm)) {
      values <- as.numeric(train_logcpm[gene, ])
      signature$training_mean[i] <- mean(values, na.rm = TRUE)
      signature$training_sd[i] <- sd(values, na.rm = TRUE)
    }
  }

  signature$usable_for_score <- signature$gene %in% rownames(train_logcpm) &
    signature$available_in_heldout &
    is.finite(signature$training_sd) & signature$training_sd > 0 &
    is.finite(signature$weight) & signature$weight > 0

  young_sig <- signature[signature$module == "young_module" & signature$usable_for_score, ]
  old_sig <- signature[signature$module == "old_module" & signature$usable_for_score, ]

  if (nrow(young_sig) == 0 || nrow(old_sig) == 0) {
    return(list(ok = FALSE, reason = "no_usable_signature_genes", signature = signature))
  }

  young_z <- (as.numeric(heldout_logcpm[young_sig$gene, heldout_sample_id]) -
    young_sig$training_mean) / young_sig$training_sd
  old_z <- (as.numeric(heldout_logcpm[old_sig$gene, heldout_sample_id]) -
    old_sig$training_mean) / old_sig$training_sd
  raw_young <- weighted_mean_safe(young_z, young_sig$weight)
  raw_old <- weighted_mean_safe(old_z, old_sig$weight)
  raw_score <- raw_young - raw_old

  train_young_scores <- rep(NA_real_, ncol(train_logcpm))
  train_old_scores <- rep(NA_real_, ncol(train_logcpm))
  train_raw_scores <- rep(NA_real_, ncol(train_logcpm))
  names(train_raw_scores) <- colnames(train_logcpm)

  for (sample_id in colnames(train_logcpm)) {
    yz <- (as.numeric(train_logcpm[young_sig$gene, sample_id]) -
      young_sig$training_mean) / young_sig$training_sd
    oz <- (as.numeric(train_logcpm[old_sig$gene, sample_id]) -
      old_sig$training_mean) / old_sig$training_sd
    train_young_scores[sample_id] <- weighted_mean_safe(yz, young_sig$weight)
    train_old_scores[sample_id] <- weighted_mean_safe(oz, old_sig$weight)
    train_raw_scores[sample_id] <- train_young_scores[sample_id] - train_old_scores[sample_id]
  }

  list(
    ok = TRUE,
    signature = signature,
    raw_young = raw_young,
    raw_old = raw_old,
    raw_score = raw_score,
    train_raw_scores = train_raw_scores
  )
}

score_train_samples <- function(signature, train_logcpm) {
  signature$training_mean <- NA_real_
  signature$training_sd <- NA_real_
  for (i in seq_len(nrow(signature))) {
    gene <- signature$gene[i]
    values <- as.numeric(train_logcpm[gene, ])
    signature$training_mean[i] <- mean(values, na.rm = TRUE)
    signature$training_sd[i] <- sd(values, na.rm = TRUE)
  }
  signature$usable_for_score <- is.finite(signature$training_sd) & signature$training_sd > 0 &
    is.finite(signature$weight) & signature$weight > 0

  young_sig <- signature[signature$module == "young_module" & signature$usable_for_score, ]
  old_sig <- signature[signature$module == "old_module" & signature$usable_for_score, ]
  raw <- rep(NA_real_, ncol(train_logcpm))
  names(raw) <- colnames(train_logcpm)
  for (sample_id in colnames(train_logcpm)) {
    yz <- (as.numeric(train_logcpm[young_sig$gene, sample_id]) -
      young_sig$training_mean) / young_sig$training_sd
    oz <- (as.numeric(train_logcpm[old_sig$gene, sample_id]) -
      old_sig$training_mean) / old_sig$training_sd
    raw[sample_id] <- weighted_mean_safe(yz, young_sig$weight) -
      weighted_mean_safe(oz, old_sig$weight)
  }
  raw
}

heldout_tmm_logcpm <- function(filtered_counts, train_metadata, heldout_metadata, train_dge) {
  train_ids <- train_metadata$sample_id
  heldout_id <- heldout_metadata$sample_id
  train_lib <- train_dge$samples$lib.size * train_dge$samples$norm.factors
  names(train_lib) <- train_ids
  ref_sample <- names(train_lib)[which.min(abs(train_lib - median(train_lib)))]
  combined_ids <- c(train_ids, heldout_id)
  combined_dge <- DGEList(counts = filtered_counts[, combined_ids, drop = FALSE])
  ref_col <- match(ref_sample, combined_ids)
  combined_dge <- calcNormFactors(combined_dge, method = "TMM", refColumn = ref_col)
  list(
    logcpm = cpm(combined_dge, log = TRUE, prior.count = 2),
    heldout_norm_factor = combined_dge$samples$norm.factors[match(heldout_id, rownames(combined_dge$samples))],
    heldout_effective_library_size = combined_dge$samples$lib.size[match(heldout_id, rownames(combined_dge$samples))] *
      combined_dge$samples$norm.factors[match(heldout_id, rownames(combined_dge$samples))],
    tmm_reference_sample = ref_sample
  )
}

build_fold_reliability <- function(outer_fit, filtered_counts, train_metadata, outer_heldout_id) {
  genes <- rownames(filtered_counts)
  de <- outer_fit$table[match(genes, outer_fit$table$gene), ]
  train_logcpm <- outer_fit$logcpm
  age_rho <- outer_fit$age_rho[genes]

  full_logfc <- de$logFC
  names(full_logfc) <- genes
  full_fdr <- de$FDR
  names(full_fdr) <- genes

  lomo_logfc <- matrix(NA_real_, nrow = length(genes), ncol = nrow(train_metadata),
                       dimnames = list(genes, train_metadata$sample_id))
  lomo_valid <- setNames(rep(FALSE, nrow(train_metadata)), train_metadata$sample_id)
  lomo_reason <- setNames(rep("", nrow(train_metadata)), train_metadata$sample_id)

  for (inner_heldout_id in train_metadata$sample_id) {
    inner_keep <- train_metadata$sample_id != inner_heldout_id
    inner_metadata <- train_metadata[inner_keep, , drop = FALSE]
    if (!is_design_full_rank(inner_metadata)) {
      lomo_reason[inner_heldout_id] <- "inner_design_not_full_rank"
      next
    }
    inner_fit <- fit_age_model(filtered_counts[, inner_metadata$sample_id, drop = FALSE], inner_metadata)
    if (!isTRUE(inner_fit$ok)) {
      lomo_reason[inner_heldout_id] <- inner_fit$reason
      next
    }
    lomo_valid[inner_heldout_id] <- TRUE
    lomo_reason[inner_heldout_id] <- "ok"
    inner_table <- inner_fit$table
    lomo_logfc[inner_table$gene, inner_heldout_id] <- inner_table$logFC
  }

  full_sign <- sign(full_logfc[genes])
  sign_match <- sweep(sign(lomo_logfc), 1, full_sign, FUN = "==")
  sign_match[full_sign == 0, ] <- TRUE
  sign_match[, !lomo_valid] <- NA
  lomo_sign_rate <- rowMeans(sign_match, na.rm = TRUE)
  valid_lomo_count <- rowSums(!is.na(sign_match))
  lomo_sign_rate[valid_lomo_count == 0] <- NA_real_

  young_inner_ids <- train_metadata$sample_id[train_metadata$age_group == "Young"]
  valid_young_inner_ids <- young_inner_ids[lomo_valid[young_inner_ids]]
  if (length(valid_young_inner_ids) > 0) {
    young_reversal <- apply(!sign_match[, valid_young_inner_ids, drop = FALSE], 1, any, na.rm = TRUE)
    young_reversal_assessed <- TRUE
  } else {
    young_reversal <- rep(FALSE, length(genes))
    young_reversal_assessed <- FALSE
  }

  low_depth_removed <- intersect(low_depth_samples, train_metadata$sample_id)
  low_depth_triggered <- length(low_depth_removed) > 0
  low_depth_valid <- FALSE
  low_depth_sign_match <- rep(NA, length(genes))
  low_depth_abs_change <- rep(NA_real_, length(genes))
  names(low_depth_sign_match) <- genes
  names(low_depth_abs_change) <- genes
  low_depth_reason <- "no_low_depth_samples_in_training"
  if (low_depth_triggered) {
    reduced_keep <- !(train_metadata$sample_id %in% low_depth_removed)
    reduced_metadata <- train_metadata[reduced_keep, , drop = FALSE]
    if (is_design_full_rank(reduced_metadata)) {
      reduced_fit <- fit_age_model(filtered_counts[, reduced_metadata$sample_id, drop = FALSE], reduced_metadata)
      if (isTRUE(reduced_fit$ok)) {
        low_depth_valid <- TRUE
        low_depth_reason <- "ok"
        reduced_logfc <- reduced_fit$table$logFC
        names(reduced_logfc) <- reduced_fit$table$gene
        low_depth_sign_match[genes] <- sign(full_logfc[genes]) == sign(reduced_logfc[genes])
        low_depth_abs_change[genes] <- abs(reduced_logfc[genes] - full_logfc[genes])
      } else {
        low_depth_reason <- reduced_fit$reason
      }
    } else {
      low_depth_reason <- "reduced_design_not_full_rank"
    }
  }
  if (!low_depth_valid) {
    low_depth_sign_match[] <- TRUE
    low_depth_abs_change[] <- NA_real_
  }

  old_male <- train_metadata$age_group == "Old" & train_metadata$sex == "male"
  old_female <- train_metadata$age_group == "Old" & train_metadata$sex == "female"
  if (sum(old_male) > 0 && sum(old_female) > 0) {
    old_male_mean <- rowMeans(train_logcpm[, old_male, drop = FALSE])
    old_female_mean <- rowMeans(train_logcpm[, old_female, drop = FALSE])
    old_sex_delta <- old_male_mean - old_female_mean
    sex_sensitive_assessed <- TRUE
  } else {
    old_sex_delta <- rep(NA_real_, length(genes))
    names(old_sex_delta) <- genes
    sex_sensitive_assessed <- FALSE
  }
  sex_effect_ratio <- abs(old_sex_delta[genes]) / pmax(abs(full_logfc[genes]), 1e-6)
  strong_sex_sensitive <- abs(old_sex_delta[genes]) > sex_sensitivity_abs_delta_threshold &
    sex_effect_ratio > sex_sensitivity_ratio_threshold
  strong_sex_sensitive[is.na(strong_sex_sensitive)] <- FALSE

  trend_compatible <- ifelse(full_logfc[genes] > 0, age_rho[genes] < 0, age_rho[genes] > 0)
  trend_compatible[is.na(trend_compatible)] <- FALSE
  sex_linked <- genes %in% sex_linked_genes

  lomo_pass <- !is.na(lomo_sign_rate) & lomo_sign_rate >= 0.9 & !young_reversal
  low_depth_pass <- as.logical(low_depth_sign_match[genes])
  low_depth_pass[is.na(low_depth_pass)] <- FALSE
  sex_pass <- !sex_linked & !strong_sex_sensitive

  reliability <- data.frame(
    gene = genes,
    adjusted_logFC = full_logfc[genes],
    FDR = full_fdr[genes],
    age_rho = age_rho[genes],
    continuous_age_trend_compatible = trend_compatible,
    LOMO_sign_rate = lomo_sign_rate,
    valid_inner_LOMO_count = valid_lomo_count,
    young_mouse_removal_direction_reversal = young_reversal,
    young_reversal_assessed = young_reversal_assessed,
    low_depth_sign_match = low_depth_pass,
    low_depth_abs_effect_size_change = low_depth_abs_change[genes],
    low_depth_triggered = low_depth_triggered,
    low_depth_valid = low_depth_valid,
    sex_effect_old_male_minus_old_female_logcpm = old_sex_delta[genes],
    sex_effect_ratio_vs_age = sex_effect_ratio,
    sex_sensitive_filtering_assessed = sex_sensitive_assessed,
    obvious_sex_linked_gene = sex_linked,
    strong_sex_sensitive_gene = strong_sex_sensitive,
    effect_size_pass_abs_logFC_gt_0_5 = abs(full_logfc[genes]) > effect_threshold,
    stringsAsFactors = FALSE
  )
  reliability$LOMO_pass <- lomo_pass
  reliability$low_depth_pass <- low_depth_pass
  reliability$sex_pass <- sex_pass
  reliability$passes_step9_initial_reliability <- reliability$effect_size_pass_abs_logFC_gt_0_5 &
    reliability$continuous_age_trend_compatible &
    reliability$LOMO_pass &
    reliability$low_depth_pass &
    reliability$sex_pass
  reliability$pi_LOMO <- pmax(0, pmin(1, reliability$LOMO_sign_rate))
  reliability$pi_depth <- ifelse(reliability$low_depth_pass, 1, 0)
  reliability$pi_sex <- ifelse(reliability$sex_pass, 1, 0)
  reliability$r_g <- reliability$pi_LOMO * reliability$pi_depth * reliability$pi_sex
  reliability$q_g <- abs(reliability$adjusted_logFC) * abs(reliability$age_rho) * reliability$r_g
  reliability$direction <- ifelse(reliability$adjusted_logFC > 0, "young_high", "old_high")

  list(
    reliability = reliability,
    lomo_valid = lomo_valid,
    lomo_reason = lomo_reason,
    low_depth_removed = low_depth_removed,
    low_depth_triggered = low_depth_triggered,
    low_depth_valid = low_depth_valid,
    low_depth_reason = low_depth_reason,
    sex_sensitive_assessed = sex_sensitive_assessed,
    young_reversal_assessed = young_reversal_assessed
  )
}

select_signature <- function(reliability, version_name, n_dir) {
  candidate <- reliability[reliability$passes_step9_initial_reliability, ]
  young <- candidate[candidate$direction == "young_high", ]
  old <- candidate[candidate$direction == "old_high", ]
  young <- young[order(-young$q_g, -abs(young$adjusted_logFC), -abs(young$age_rho), young$gene), ]
  old <- old[order(-old$q_g, -abs(old$adjusted_logFC), -abs(old$age_rho), old$gene), ]
  selected <- rbind(head(young, n_dir), head(old, n_dir))
  selected$version <- version_name
  selected$genes_per_direction_target <- n_dir
  selected$module <- ifelse(selected$direction == "young_high", "young_module", "old_module")
  selected$weight <- pmin(abs(selected$adjusted_logFC), w_max) * selected$r_g
  selected
}

jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  if (length(union(a, b)) == 0) {
    return(NA_real_)
  }
  length(intersect(a, b)) / length(union(a, b))
}

message("Reading inputs")
counts_raw <- readRDS(counts_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(colnames(counts_raw), metadata$sample_id)) {
  stop("Raw pseudobulk count columns do not match metadata sample_id order")
}

score_rows <- list()
signature_rows <- list()
diagnostic_rows <- list()

for (heldout_id in metadata$sample_id) {
  message(sprintf("Outer fold held out: %s", heldout_id))
  train_keep <- metadata$sample_id != heldout_id
  train_metadata <- metadata[train_keep, , drop = FALSE]
  heldout_metadata <- metadata[!train_keep, , drop = FALSE]
  fold_type <- ifelse(heldout_metadata$age_group == "Young", "held_out_young_stress_test", "held_out_old")

  design <- design_for(train_metadata)
  design_full_rank <- qr(design)$rank == ncol(design)
  gene_keep <- filter_genes_train(counts_raw[, train_keep, drop = FALSE])
  filtered_counts <- counts_raw[gene_keep, , drop = FALSE]

  outer_fit <- NULL
  reliability_result <- NULL
  if (design_full_rank) {
    outer_fit <- fit_age_model(filtered_counts[, train_keep, drop = FALSE], train_metadata)
  } else {
    outer_fit <- list(ok = FALSE, reason = "outer_design_not_full_rank")
  }

  if (isTRUE(outer_fit$ok)) {
    reliability_result <- build_fold_reliability(outer_fit, filtered_counts, train_metadata, heldout_id)
    heldout_norm <- heldout_tmm_logcpm(filtered_counts, train_metadata, heldout_metadata, outer_fit$dge)
  }

  for (v in seq_len(nrow(versions))) {
    version_name <- versions$version[v]
    n_dir <- versions$genes_per_direction[v]

    if (!isTRUE(outer_fit$ok)) {
      diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
        held_out_mouse = heldout_id,
        version = version_name,
        fold_type = fold_type,
        design_full_rank = design_full_rank,
        fold_ok = FALSE,
        failure_reason = outer_fit$reason,
        n_training_mice = nrow(train_metadata),
        n_training_young = sum(train_metadata$age_group == "Young"),
        n_training_old = sum(train_metadata$age_group == "Old"),
        genes_after_fold_filter = sum(gene_keep),
        candidate_genes_total = NA_integer_,
        candidate_young_high = NA_integer_,
        candidate_old_high = NA_integer_,
        selected_young_high = NA_integer_,
        selected_old_high = NA_integer_,
        low_depth_triggered = NA,
        low_depth_removed = NA_character_,
        low_depth_valid = NA,
        low_depth_reason = NA_character_,
        sex_sensitive_filtering_assessed = NA,
        young_reversal_assessed = NA,
        valid_inner_lomo_folds = NA_integer_,
        invalid_inner_lomo_folds = NA_integer_,
        tmm_reference_sample = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }

    reliability <- reliability_result$reliability
    signature <- select_signature(reliability, version_name, n_dir)
    signature$held_out_mouse <- heldout_id
    signature$fold_type <- fold_type
    signature_rows[[length(signature_rows) + 1]] <- signature

    candidate <- reliability[reliability$passes_step9_initial_reliability, ]
    candidate_young <- sum(candidate$direction == "young_high")
    candidate_old <- sum(candidate$direction == "old_high")
    selected_young <- sum(signature$direction == "young_high")
    selected_old <- sum(signature$direction == "old_high")

    train_raw <- score_train_samples(signature, outer_fit$logcpm)
    train_age_group <- train_metadata$age_group[match(names(train_raw), train_metadata$sample_id)]
    young_center <- median(train_raw[train_age_group == "Young"], na.rm = TRUE)
    old_center <- median(train_raw[train_age_group == "Old"], na.rm = TRUE)
    calibration_denominator <- young_center - old_center

    score_result <- score_one_sample(signature, outer_fit$logcpm, heldout_norm$logcpm, heldout_id)
    if (isTRUE(score_result$ok) && is.finite(calibration_denominator) &&
        abs(calibration_denominator) > .Machine$double.eps) {
      predicted_raw_score <- score_result$raw_score
      predicted_young_raw <- score_result$raw_young
      predicted_old_raw <- score_result$raw_old
      predicted_calibrated <- (score_result$raw_score - old_center) / calibration_denominator
      total_weight <- sum(abs(signature$weight[score_result$signature$usable_for_score]))
      available <- score_result$signature$usable_for_score & score_result$signature$available_in_heldout
      young_available <- available & score_result$signature$module == "young_module"
      old_available <- available & score_result$signature$module == "old_module"
      heldout_gene_coverage <- sum(available) / nrow(signature)
      heldout_weighted_coverage <- sum(abs(score_result$signature$weight[available])) / total_weight
      heldout_young_coverage <- sum(young_available) / sum(score_result$signature$module == "young_module")
      heldout_old_coverage <- sum(old_available) / sum(score_result$signature$module == "old_module")
      score_ok <- TRUE
      score_reason <- "ok"
    } else {
      predicted_raw_score <- NA_real_
      predicted_young_raw <- NA_real_
      predicted_old_raw <- NA_real_
      predicted_calibrated <- NA_real_
      heldout_gene_coverage <- NA_real_
      heldout_weighted_coverage <- NA_real_
      heldout_young_coverage <- NA_real_
      heldout_old_coverage <- NA_real_
      score_ok <- FALSE
      score_reason <- if (isTRUE(score_result$ok)) {
        "invalid_calibration_denominator"
      } else if (!is.null(score_result$reason)) {
        score_result$reason
      } else {
        "score_failed"
      }
    }

    score_rows[[length(score_rows) + 1]] <- data.frame(
      mouse = heldout_id,
      age = heldout_metadata$age_months,
      sex = as.character(heldout_metadata$sex),
      version = version_name,
      fold_type = fold_type,
      predicted_raw_score = predicted_raw_score,
      predicted_young_module_raw = predicted_young_raw,
      predicted_old_module_raw = predicted_old_raw,
      predicted_calibrated_score = predicted_calibrated,
      signature_size = nrow(signature),
      selected_young_high = selected_young,
      selected_old_high = selected_old,
      design_full_rank = design_full_rank,
      score_ok = score_ok,
      score_failure_reason = score_reason,
      heldout_gene_coverage = heldout_gene_coverage,
      heldout_young_module_coverage = heldout_young_coverage,
      heldout_old_module_coverage = heldout_old_coverage,
      heldout_weighted_coverage = heldout_weighted_coverage,
      calibration_denominator = calibration_denominator,
      young_reference_center = young_center,
      old_reference_center = old_center,
      fold_candidate_genes_total = nrow(candidate),
      fold_candidate_young_high = candidate_young,
      fold_candidate_old_high = candidate_old,
      low_depth_triggered = reliability_result$low_depth_triggered,
      low_depth_removed = paste(reliability_result$low_depth_removed, collapse = ";"),
      low_depth_valid = reliability_result$low_depth_valid,
      sex_sensitive_filtering = reliability_result$sex_sensitive_assessed,
      heldout_tmm_norm_factor = heldout_norm$heldout_norm_factor,
      heldout_effective_library_size = heldout_norm$heldout_effective_library_size,
      heldout_cell_count = heldout_metadata$cell_count,
      stringsAsFactors = FALSE
    )

    diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
      held_out_mouse = heldout_id,
      version = version_name,
      fold_type = fold_type,
      design_full_rank = design_full_rank,
      fold_ok = TRUE,
      failure_reason = "ok",
      n_training_mice = nrow(train_metadata),
      n_training_young = sum(train_metadata$age_group == "Young"),
      n_training_old = sum(train_metadata$age_group == "Old"),
      genes_after_fold_filter = sum(gene_keep),
      candidate_genes_total = nrow(candidate),
      candidate_young_high = candidate_young,
      candidate_old_high = candidate_old,
      selected_young_high = selected_young,
      selected_old_high = selected_old,
      low_depth_triggered = reliability_result$low_depth_triggered,
      low_depth_removed = paste(reliability_result$low_depth_removed, collapse = ";"),
      low_depth_valid = reliability_result$low_depth_valid,
      low_depth_reason = reliability_result$low_depth_reason,
      sex_sensitive_filtering_assessed = reliability_result$sex_sensitive_assessed,
      young_reversal_assessed = reliability_result$young_reversal_assessed,
      valid_inner_lomo_folds = sum(reliability_result$lomo_valid),
      invalid_inner_lomo_folds = sum(!reliability_result$lomo_valid),
      tmm_reference_sample = heldout_norm$tmm_reference_sample,
      stringsAsFactors = FALSE
    )
  }
}

scores <- do.call(rbind, score_rows)
fold_signatures <- do.call(rbind, signature_rows)
diagnostics <- do.call(rbind, diagnostic_rows)

overlap_rows <- list()
for (version_name in versions$version) {
  version_sigs <- fold_signatures[fold_signatures$version == version_name, ]
  fold_ids <- unique(version_sigs$held_out_mouse)
  for (i in seq_along(fold_ids)) {
    for (j in seq_along(fold_ids)) {
      if (j <= i) next
      a <- version_sigs[version_sigs$held_out_mouse == fold_ids[i], ]
      b <- version_sigs[version_sigs$held_out_mouse == fold_ids[j], ]
      overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
        version = version_name,
        fold_a = fold_ids[i],
        fold_b = fold_ids[j],
        jaccard_all = jaccard(a$gene, b$gene),
        jaccard_young_module = jaccard(a$gene[a$module == "young_module"], b$gene[b$module == "young_module"]),
        jaccard_old_module = jaccard(a$gene[a$module == "old_module"], b$gene[b$module == "old_module"]),
        size_a = nrow(a),
        size_b = nrow(b),
        stringsAsFactors = FALSE
      )
    }
  }
}
fold_overlap <- do.call(rbind, overlap_rows)

summary_rows <- list()
for (version_name in versions$version) {
  x <- scores[scores$version == version_name & scores$score_ok, ]
  old_x <- x[x$fold_type == "held_out_old", ]
  old_x$sex_numeric <- ifelse(old_x$sex == "male", 1, 0)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    version = version_name,
    n_scored_folds = nrow(x),
    n_old_folds = nrow(old_x),
    n_young_stress_folds = sum(x$fold_type == "held_out_young_stress_test"),
    spearman_score_vs_age_all = suppressWarnings(cor(x$predicted_calibrated_score, x$age, method = "spearman")),
    spearman_score_vs_age_old_only = suppressWarnings(cor(old_x$predicted_calibrated_score, old_x$age, method = "spearman")),
    spearman_score_vs_effective_library_size_all = suppressWarnings(cor(x$predicted_calibrated_score, x$heldout_effective_library_size, method = "spearman")),
    spearman_score_vs_cell_count_all = suppressWarnings(cor(x$predicted_calibrated_score, x$heldout_cell_count, method = "spearman")),
    old_male_median_score = median(old_x$predicted_calibrated_score[old_x$sex == "male"], na.rm = TRUE),
    old_female_median_score = median(old_x$predicted_calibrated_score[old_x$sex == "female"], na.rm = TRUE),
    old_male_minus_female_median_score = median(old_x$predicted_calibrated_score[old_x$sex == "male"], na.rm = TRUE) -
      median(old_x$predicted_calibrated_score[old_x$sex == "female"], na.rm = TRUE),
    young_stress_min_score = min(x$predicted_calibrated_score[x$fold_type == "held_out_young_stress_test"], na.rm = TRUE),
    old_fold_median_score = median(old_x$predicted_calibrated_score, na.rm = TRUE),
    old_folds_below_young_stress_min = sum(old_x$predicted_calibrated_score < min(x$predicted_calibrated_score[x$fold_type == "held_out_young_stress_test"], na.rm = TRUE), na.rm = TRUE),
    median_signature_jaccard_all = median(fold_overlap$jaccard_all[fold_overlap$version == version_name], na.rm = TRUE),
    median_signature_jaccard_young = median(fold_overlap$jaccard_young_module[fold_overlap$version == version_name], na.rm = TRUE),
    median_signature_jaccard_old = median(fold_overlap$jaccard_old_module[fold_overlap$version == version_name], na.rm = TRUE),
    min_heldout_gene_coverage = min(x$heldout_gene_coverage, na.rm = TRUE),
    min_heldout_weighted_coverage = min(x$heldout_weighted_coverage, na.rm = TRUE),
    interpretation = "nested LOMO; young-heldout folds are stress tests because training has one young mouse",
    stringsAsFactors = FALSE
  )
}
summary_table <- do.call(rbind, summary_rows)

write.csv(scores, file.path(validation_dir, "step14_nested_lomo_scores.csv"), row.names = FALSE)
write.csv(fold_signatures, file.path(validation_dir, "step14_fold_signatures.csv"), row.names = FALSE)
write.csv(fold_overlap, file.path(validation_dir, "step14_fold_signature_overlap.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(validation_dir, "step14_fold_diagnostics.csv"), row.names = FALSE)
write.csv(summary_table, file.path(validation_dir, "step14_nested_lomo_summary.csv"), row.names = FALSE)

report_lines <- c(
  "# Step 14: Nested Leave-One-Mouse-Out Score Validation",
  "",
  "## Core Leakage Control",
  "",
  "Each outer fold removed one mouse before any fold-specific normalization, DE, stability filtering, ranking, signature selection, gene standardization, or score calibration. The full-data Step 10/12 signatures were not reused for held-out scoring.",
  "",
  "## Per-Fold Training Workflow",
  "",
  "1. Removed one mouse.",
  "2. Re-applied the CPM gene filter on the remaining training mice.",
  "3. Recomputed TMM normalization on the training mice.",
  "4. Refit the sex-adjusted edgeR model `~ sex + age_group`.",
  "5. Recomputed continuous age Spearman correlations on training logCPM.",
  "6. Recomputed inner training-side LOMO stability where the design remained full rank.",
  "7. Recomputed the prespecified low-depth sensitivity by removing available low-depth training samples.",
  "8. Recomputed sex-linked and old-sex-sensitive filtering.",
  "9. Re-ranked candidates with `q_g = |adjusted_logFC| * |age_rho| * r_g`.",
  "10. Re-selected Small/Medium/Large signatures inside the fold.",
  "11. Recomputed `mu_g`, `sd_g`, `M_Y`, and `M_O` on training mice only.",
  "12. Normalized the held-out mouse using a TMM factor against a training reference sample, then scored it.",
  "",
  "## Known Limitation",
  "",
  "When a young mouse is held out, the training set contains only one young mouse. Those folds are marked `held_out_young_stress_test`; they can run but should not be treated as stable generalization estimates.",
  "",
  "## Summary",
  "",
  paste(capture.output(print(summary_table, row.names = FALSE)), collapse = "\n"),
  "",
  "## Design Diagnostics",
  "",
  sprintf("- Outer folds with full-rank `~ sex + age_group` design: %s / %s",
          sum(unique(diagnostics[, c("held_out_mouse", "design_full_rank")])$design_full_rank),
          length(unique(diagnostics$held_out_mouse))),
  sprintf("- Fold-version rows with successful scoring: %s / %s", sum(scores$score_ok), nrow(scores)),
  sprintf("- Minimum held-out gene coverage: %.3f", min(scores$heldout_gene_coverage, na.rm = TRUE)),
  sprintf("- Minimum held-out weighted coverage: %.3f", min(scores$heldout_weighted_coverage, na.rm = TRUE)),
  "",
  "## Outputs",
  "",
  "- `outputs/validation/step14_nested_lomo_scores.csv`",
  "- `outputs/validation/step14_fold_signatures.csv`",
  "- `outputs/validation/step14_fold_signature_overlap.csv`",
  "- `outputs/validation/step14_fold_diagnostics.csv`",
  "- `outputs/validation/step14_nested_lomo_summary.csv`"
)
writeLines(report_lines, file.path(validation_dir, "step14_nested_lomo_report.md"))

cat("Step 14 nested LOMO complete\n")
print(summary_table, row.names = FALSE)
