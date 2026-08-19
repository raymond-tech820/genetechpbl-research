#!/usr/bin/env Rscript

set.seed(20260721)

root <- getwd()
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

out_root <- file.path(root, "outputs", "facs_v2")
processed_dir <- file.path(out_root, "processed")
out_dir <- file.path(out_root, "technical_burden")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.rds"))
metadata <- read.csv(file.path(processed_dir, "facs_v2_limb_msc_mouse_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(colnames(counts), metadata$mouse))

metadata$age_group <- factor(metadata$age_group, levels = c("Young", "Old"))
metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
metadata$age_sex_group <- interaction(metadata$sex, metadata$age_group, sep = "_", drop = TRUE)

cpm_threshold <- 1
min_mice_expressed <- 2
signature_sizes <- c(Small = 20, Medium = 50, Large = 100)
lambda_grid <- c(original = 0, tech_mild = 1, tech_moderate = 2, tech_strong = 4)
epsilon <- 1e-6
sex_linked_genes <- c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Rps4y1", "Rps4y2", "Sry", "Jarid1d")
cell_cycle_genes <- unique(c("Mki67", "Top2a", "Pcna", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6", "Mcm7", "Ccna2", "Ccnb1", "Ccnb2", "Ccnd1", "Ccne1", "Cdk1", "Cdk2", "Cdk4", "Birc5", "Ube2c", "Cenpf", "Nusap1", "Tyms", "Tk1", "Hmgb2"))

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

auc_rank <- function(labels, scores) {
  labels <- as.integer(labels)
  ok <- is.finite(scores) & !is.na(labels)
  labels <- labels[ok]; scores <- scores[ok]
  n_pos <- sum(labels == 1); n_neg <- sum(labels == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

compute_residual_burden <- function(logcpm, design, train_meta) {
  fit <- lmFit(logcpm, design)
  fitted <- fit$coefficients %*% t(design)
  residuals <- logcpm - fitted
  raw_lib <- log1p(train_meta$pseudobulk_library_size)
  cell_count <- log1p(train_meta$n_cells)
  detected <- train_meta$pseudobulk_detected_genes
  eff_lib <- log1p(train_meta$effective_library_size_tmm)
  burden <- t(vapply(rownames(residuals), function(g) {
    e <- as.numeric(residuals[g, ])
    x <- as.numeric(logcpm[g, ])
    c(
      residual_cor_library = abs(safe_cor(e, raw_lib)),
      residual_cor_cell_count = abs(safe_cor(e, cell_count)),
      residual_cor_detected = abs(safe_cor(e, detected)),
      residual_cor_effective_library = abs(safe_cor(e, eff_lib)),
      marginal_cor_library = abs(safe_cor(x, raw_lib)),
      marginal_cor_cell_count = abs(safe_cor(x, cell_count)),
      marginal_cor_detected = abs(safe_cor(x, detected))
    )
  }, numeric(7)))
  burden <- as.data.frame(burden)
  burden$gene <- rownames(residuals)
  burden$technical_burden_mean <- rowMeans(burden[, c("residual_cor_library", "residual_cor_cell_count", "residual_cor_detected")], na.rm = TRUE)
  burden$technical_burden_max <- apply(burden[, c("residual_cor_library", "residual_cor_cell_count", "residual_cor_detected")], 1, max, na.rm = TRUE)
  burden
}

build_signature <- function(ranked, branch, per_direction) {
  ranked$branch_rank_score <- ranked[[paste0("rank_score_", branch)]]
  young <- ranked[ranked$module == "young_high" & ranked$reliability_pass, , drop = FALSE]
  old <- ranked[ranked$module == "old_high" & ranked$reliability_pass, , drop = FALSE]
  young <- young[order(-young$branch_rank_score, young$gene), , drop = FALSE]
  old <- old[order(-old$branch_rank_score, old$gene), , drop = FALSE]
  sig <- rbind(head(young, per_direction), head(old, per_direction))
  sig$branch <- branch
  sig$signature_target_per_direction <- per_direction
  sig
}

score_signature <- function(sig, train_expr, held_expr, train_meta) {
  genes <- intersect(sig$gene, rownames(train_expr))
  genes <- intersect(genes, rownames(held_expr))
  sig <- sig[match(genes, sig$gene), , drop = FALSE]
  if (length(genes) == 0 || nrow(sig) == 0) return(rep(NA_real_, 7))
  mu <- rowMeans(train_expr[genes, , drop = FALSE])
  sdv <- apply(train_expr[genes, , drop = FALSE], 1, sd)
  usable <- is.finite(mu) & is.finite(sdv) & sdv > 0
  genes <- genes[usable]; sig <- sig[usable, , drop = FALSE]; mu <- mu[usable]; sdv <- sdv[usable]
  if (length(genes) == 0) return(rep(NA_real_, 7))
  z_train <- sweep(sweep(train_expr[genes, , drop = FALSE], 1, mu, "-"), 1, sdv, "/")
  z_held <- (as.numeric(held_expr[genes, 1]) - mu) / sdv
  w <- sig$weight
  w[!is.finite(w) | w <= 0] <- 1
  module_score <- function(z, module) {
    idx <- which(sig$module == module)
    if (length(idx) == 0) return(NA_real_)
    weighted.mean(z[idx], w[idx], na.rm = TRUE)
  }
  train_raw <- vapply(seq_len(ncol(z_train)), function(j) {
    module_score(z_train[, j], "young_high") - module_score(z_train[, j], "old_high")
  }, numeric(1))
  held_raw <- module_score(z_held, "young_high") - module_score(z_held, "old_high")
  young_center <- median(train_raw[train_meta$age_group == "Young"], na.rm = TRUE)
  old_center <- median(train_raw[train_meta$age_group == "Old"], na.rm = TRUE)
  denom <- young_center - old_center
  relative <- if (is.finite(denom) && abs(denom) > epsilon) (held_raw - old_center) / denom else NA_real_
  total_w <- sum(sig$weight, na.rm = TRUE)
  c(raw = held_raw, relative = relative, coverage = nrow(sig) / nrow(sig), weighted_coverage = 1,
    young_center = young_center, old_center = old_center, denom = denom)
}

run_nested <- function(meta, count_mat, analysis_label) {
  score_rows <- list(); sig_rows <- list(); gene_rows <- list(); diag_rows <- list()
  for (heldout in meta$mouse) {
    message("[", analysis_label, "] held out: ", heldout)
    held_idx <- which(meta$mouse == heldout)
    train_idx <- setdiff(seq_len(nrow(meta)), held_idx)
    train_meta <- meta[train_idx, , drop = FALSE]
    held_meta <- meta[held_idx, , drop = FALSE]
    train_meta$age_group <- factor(train_meta$age_group, levels = c("Young", "Old"))
    train_meta$sex <- factor(train_meta$sex, levels = c("female", "male"))
    train_meta$age_sex_group <- factor(interaction(train_meta$sex, train_meta$age_group, sep = "_", drop = TRUE),
      levels = c("female_Young", "male_Young", "female_Old", "male_Old"))
    support <- table(train_meta$age_sex_group)
    weak_support <- any(support[c("female_Young", "male_Young", "female_Old", "male_Old")] < 2)

    train_raw <- count_mat[, train_meta$mouse, drop = FALSE]
    held_raw <- count_mat[, held_meta$mouse, drop = FALSE]
    keep <- rowSums(t(t(train_raw) / colSums(train_raw)) * 1e6 > cpm_threshold) >= min_mice_expressed
    train_counts <- train_raw[keep, , drop = FALSE]
    held_counts <- held_raw[keep, , drop = FALSE]

    dge <- DGEList(counts = train_counts)
    dge <- calcNormFactors(dge, method = "TMM")
    eff <- dge$samples$lib.size * dge$samples$norm.factors
    train_meta$effective_library_size_tmm <- as.numeric(eff)

    design_age <- model.matrix(~ age_group, data = train_meta)
    design_fact <- model.matrix(~ 0 + age_sex_group, data = train_meta)
    full_rank <- qr(design_fact)$rank == ncol(design_fact)

    v_age <- voom(dge, design_age, plot = FALSE)
    fit_age <- eBayes(lmFit(v_age, design_age))
    age_tab <- topTable(fit_age, coef = "age_groupOld", number = Inf, sort.by = "none")

    v_fact <- voom(dge, design_fact, plot = FALSE)
    fit_fact <- lmFit(v_fact, design_fact)
    cm <- makeContrasts(
      female_age = age_sex_groupfemale_Old - age_sex_groupfemale_Young,
      male_age = age_sex_groupmale_Old - age_sex_groupmale_Young,
      common_age = ((age_sex_groupfemale_Old - age_sex_groupfemale_Young) + (age_sex_groupmale_Old - age_sex_groupmale_Young)) / 2,
      interaction_age_by_sex = (age_sex_groupmale_Old - age_sex_groupmale_Young) - (age_sex_groupfemale_Old - age_sex_groupfemale_Young),
      levels = design_fact
    )
    fit_con <- eBayes(contrasts.fit(fit_fact, cm))
    female_tab <- topTable(fit_con, coef = "female_age", number = Inf, sort.by = "none")
    male_tab <- topTable(fit_con, coef = "male_age", number = Inf, sort.by = "none")
    common_tab <- topTable(fit_con, coef = "common_age", number = Inf, sort.by = "none")
    interaction_tab <- topTable(fit_con, coef = "interaction_age_by_sex", number = Inf, sort.by = "none")

    genes <- rownames(train_counts)
    rownames(age_tab) <- rownames(female_tab) <- rownames(male_tab) <- rownames(common_tab) <- rownames(interaction_tab) <- genes
    ranked <- data.frame(
      gene = genes,
      age_only_logFC = age_tab[genes, "logFC"],
      age_only_t = age_tab[genes, "t"],
      female_logFC = female_tab[genes, "logFC"],
      female_t = female_tab[genes, "t"],
      male_logFC = male_tab[genes, "logFC"],
      male_t = male_tab[genes, "t"],
      common_logFC = common_tab[genes, "logFC"],
      common_t = common_tab[genes, "t"],
      interaction_logFC = interaction_tab[genes, "logFC"],
      stringsAsFactors = FALSE
    )
    ranked$sex_linked_flag <- ranked$gene %in% sex_linked_genes
    ranked$female_direction <- sign(ranked$female_logFC)
    ranked$male_direction <- sign(ranked$male_logFC)
    ranked$common_direction <- sign(ranked$common_logFC)
    ranked$age_only_direction <- sign(ranked$age_only_logFC)
    ranked$sex_direction_concordant <- ranked$female_direction == ranked$male_direction & ranked$female_direction != 0
    ranked$age_only_common_concordant <- ranked$age_only_direction == ranked$common_direction & ranked$common_direction != 0
    ranked$interaction_magnitude_ratio <- abs(ranked$interaction_logFC) / (abs(ranked$female_logFC) + abs(ranked$male_logFC) + epsilon)
    ranked$sex_consistency_penalty <- 1 / (1 + ranked$interaction_magnitude_ratio)
    ranked$module <- ifelse(ranked$common_logFC < 0, "young_high", ifelse(ranked$common_logFC > 0, "old_high", "neutral"))
    ranked$reliability_pass <- with(ranked, sex_direction_concordant & age_only_common_concordant & !sex_linked_flag & is.finite(common_logFC) & is.finite(common_t))
    ranked$weight <- abs(ranked$common_logFC) * abs(ranked$common_t) * ranked$sex_consistency_penalty
    ranked$rank_score_original <- ranked$weight * as.numeric(ranked$reliability_pass)

    logcpm <- cpm(dge, log = TRUE, prior.count = 1)
    burden <- compute_residual_burden(logcpm, design_fact, train_meta)
    ranked <- merge(ranked, burden, by = "gene", sort = FALSE)
    for (branch in names(lambda_grid)) {
      lambda <- lambda_grid[[branch]]
      ranked[[paste0("technical_penalty_", branch)]] <- 1 / (1 + lambda * ranked$technical_burden_mean)
      ranked[[paste0("rank_score_", branch)]] <- ranked$rank_score_original * ranked[[paste0("technical_penalty_", branch)]]
    }

    deploy_train <- log2(t(t(train_counts) / colSums(train_counts)) * 1e6 + 1)
    deploy_held <- log2(t(t(held_counts) / colSums(held_counts)) * 1e6 + 1)

    gene_rows[[paste(analysis_label, heldout)]] <- transform(
      ranked,
      analysis = analysis_label,
      heldout_mouse = heldout,
      weak_support_fold = weak_support
    )
    diag_rows[[paste(analysis_label, heldout)]] <- data.frame(
      analysis = analysis_label,
      heldout_mouse = heldout,
      heldout_age = held_meta$age,
      heldout_age_group = as.character(held_meta$age_group),
      heldout_sex = as.character(held_meta$sex),
      design_full_rank = full_rank,
      weak_support_fold = weak_support,
      retained_genes = sum(keep),
      reliability_pass_genes = sum(ranked$reliability_pass),
      stringsAsFactors = FALSE
    )

    for (branch in names(lambda_grid)) {
      for (version in names(signature_sizes)) {
        sig <- build_signature(ranked, branch, signature_sizes[[version]])
        sig$analysis <- analysis_label
        sig$heldout_mouse <- heldout
        sig$signature_version <- version
        sig$mitochondrial_flag <- grepl("^mt-|^Mt-", sig$gene)
        sig$ribosomal_flag <- grepl("^Rp[sl]", sig$gene)
        sig$cell_cycle_flag <- sig$gene %in% cell_cycle_genes
        sig_rows[[paste(analysis_label, heldout, branch, version)]] <- sig
        sc <- score_signature(sig, deploy_train, deploy_held, train_meta)
        score_rows[[paste(analysis_label, heldout, branch, version)]] <- data.frame(
          analysis = analysis_label,
          mouse = heldout,
          age = held_meta$age,
          age_months = held_meta$age_months,
          age_group = as.character(held_meta$age_group),
          sex = as.character(held_meta$sex),
          weak_support_fold = weak_support,
          ranking_branch = branch,
          lambda = lambda_grid[[branch]],
          signature_version = version,
          signature_size = nrow(sig),
          young_high_n = sum(sig$module == "young_high"),
          old_high_n = sum(sig$module == "old_high"),
          predicted_raw_score = sc[["raw"]],
          predicted_calibrated_score = sc[["relative"]],
          gene_coverage = sc[["coverage"]],
          weighted_gene_coverage = sc[["weighted_coverage"]],
          young_reference_center = sc[["young_center"]],
          old_reference_center = sc[["old_center"]],
          calibration_denominator = sc[["denom"]],
          median_technical_burden = median(sig$technical_burden_mean, na.rm = TRUE),
          mean_technical_burden = mean(sig$technical_burden_mean, na.rm = TRUE),
          mitochondrial_fraction = mean(sig$mitochondrial_flag),
          ribosomal_fraction = mean(sig$ribosomal_flag),
          cell_cycle_fraction = mean(sig$cell_cycle_flag),
          heldout_raw_library_size = held_meta$pseudobulk_library_size,
          heldout_detected_genes = held_meta$pseudobulk_detected_genes,
          heldout_cell_count = held_meta$n_cells,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    scores = do.call(rbind, score_rows),
    signatures = do.call(rbind, sig_rows),
    genes = do.call(rbind, gene_rows),
    diagnostics = do.call(rbind, diag_rows)
  )
}

summarize_scores <- function(scores, signatures) {
  rows <- list()
  for (analysis in unique(scores$analysis)) {
    for (branch in unique(scores$ranking_branch)) {
      for (version in unique(scores$signature_version)) {
        sv <- scores[scores$analysis == analysis & scores$ranking_branch == branch & scores$signature_version == version, ]
        sigv <- signatures[signatures$analysis == analysis & signatures$branch == branch & signatures$signature_version == version, ]
        sets <- split(sigv$gene, sigv$heldout_mouse)
        jac <- NA_real_
        if (length(sets) > 1) {
          pairs <- combn(names(sets), 2, simplify = FALSE)
          jac <- median(vapply(pairs, function(p) length(intersect(sets[[p[1]]], sets[[p[2]]])) / length(union(sets[[p[1]]], sets[[p[2]]])), numeric(1)), na.rm = TRUE)
        }
        freq <- aggregate(heldout_mouse ~ gene + module, sigv, length)
        genes_75 <- sum(freq$heldout_mouse / length(unique(sv$mouse)) >= 0.75)
        rows[[paste(analysis, branch, version)]] <- data.frame(
          analysis = analysis,
          ranking_branch = branch,
          lambda = unique(sv$lambda),
          signature_version = version,
          n_folds = nrow(sv),
          zero_signature_folds = sum(sv$signature_size == 0),
          auc = auc_rank(as.integer(sv$age_group == "Young"), sv$predicted_calibrated_score),
          all_age_rho = safe_cor(sv$predicted_calibrated_score, sv$age_months),
          old_only_rho = safe_cor(sv$predicted_calibrated_score[sv$age_group == "Old"], sv$age_months[sv$age_group == "Old"]),
          young_minus_old_median = median(sv$predicted_calibrated_score[sv$age_group == "Young"], na.rm = TRUE) - median(sv$predicted_calibrated_score[sv$age_group == "Old"], na.rm = TRUE),
          library_rho = safe_cor(sv$predicted_calibrated_score, sv$heldout_raw_library_size),
          detected_rho = safe_cor(sv$predicted_calibrated_score, sv$heldout_detected_genes),
          cell_count_rho = safe_cor(sv$predicted_calibrated_score, sv$heldout_cell_count),
          median_jaccard = jac,
          genes_selected_ge_75pct = genes_75,
          median_signature_technical_burden = median(sv$median_technical_burden, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

message("Running full nested technical-burden sensitivity")
main <- run_nested(metadata, counts, "all_mice")
lowest_depth_mouse <- metadata$mouse[which.min(metadata$pseudobulk_library_size)]
metadata_no_low <- metadata[metadata$mouse != lowest_depth_mouse, , drop = FALSE]
counts_no_low <- counts[, metadata_no_low$mouse, drop = FALSE]
no_low <- run_nested(metadata_no_low, counts_no_low, paste0("remove_lowest_depth_", lowest_depth_mouse))

scores <- rbind(main$scores, no_low$scores)
signatures <- rbind(main$signatures, no_low$signatures)
genes <- rbind(main$genes, no_low$genes)
diagnostics <- rbind(main$diagnostics, no_low$diagnostics)
summary <- summarize_scores(scores, signatures)

write.csv(genes, file.path(out_dir, "technical_burden_gene_metrics_by_fold.csv"), row.names = FALSE)
write.csv(signatures, file.path(out_dir, "technical_penalty_fold_signatures.csv"), row.names = FALSE)
write.csv(scores, file.path(out_dir, "technical_penalty_nested_lomo_scores.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out_dir, "technical_penalty_fold_diagnostics.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "technical_penalty_model_summary.csv"), row.names = FALSE)

sel_freq <- aggregate(heldout_mouse ~ analysis + branch + signature_version + gene + module, signatures, length)
names(sel_freq)[names(sel_freq) == "heldout_mouse"] <- "selected_folds"
fold_counts <- aggregate(mouse ~ analysis + ranking_branch + signature_version, scores, length)
names(fold_counts) <- c("analysis", "branch", "signature_version", "n_folds")
sel_freq <- merge(sel_freq, fold_counts, by = c("analysis", "branch", "signature_version"), all.x = TRUE)
sel_freq$selection_frequency <- sel_freq$selected_folds / sel_freq$n_folds
write.csv(sel_freq[order(sel_freq$analysis, sel_freq$branch, sel_freq$signature_version, -sel_freq$selection_frequency, sel_freq$gene), ],
          file.path(out_dir, "technical_penalty_gene_selection_frequency.csv"), row.names = FALSE)

report <- c(
  "# FACS v2 Step 12: Technical-Burden Sensitivity Audit",
  "",
  "## Scope",
  "",
  "The original factorial ranking is retained as the primary specification. Technical-burden penalties are sensitivity branches only.",
  "",
  "Technical burden is estimated within each outer training fold from residual expression after `~0 + sex:age_group` adjustment. The held-out mouse is not used to estimate burden, ranking, weights, or calibration.",
  "",
  "No files under `data_facs` were modified.",
  "",
  "## Penalty",
  "",
  "Composite burden: mean absolute residual Spearman correlation with raw library size, cell count, and detected genes.",
  "",
  "Penalty: `1 / (1 + lambda * burden)`.",
  "",
  "Branches: original lambda=0, tech_mild lambda=1, tech_moderate lambda=2, tech_strong lambda=4.",
  "",
  paste0("Lowest-depth donor full rerun excludes: `", lowest_depth_mouse, "`."),
  "",
  "## Summary",
  "",
  paste(capture.output(print(summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Outputs",
  "",
  "- `technical_burden_gene_metrics_by_fold.csv`",
  "- `technical_penalty_fold_signatures.csv`",
  "- `technical_penalty_nested_lomo_scores.csv`",
  "- `technical_penalty_fold_diagnostics.csv`",
  "- `technical_penalty_model_summary.csv`",
  "- `technical_penalty_gene_selection_frequency.csv`"
)
writeLines(report, file.path(out_dir, "technical_burden_sensitivity_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step12.txt"))

message("Done")
message(sprintf("Report: %s", file.path(out_dir, "technical_burden_sensitivity_report.md")))
