#!/usr/bin/env Rscript

set.seed(20260721)

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
validation_dir <- file.path(out_root, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.rds")
metadata_path <- file.path(processed_dir, "facs_v2_limb_msc_mouse_metadata.csv")

cpm_threshold <- 1
min_mice_expressed <- 2
signature_sizes <- c(Small = 20, Medium = 50, Large = 100)
epsilon <- 1e-6
top_variable_n <- 2000

sex_linked_genes <- c(
  "Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2",
  "Rps4y1", "Rps4y2", "Sry", "Jarid1d"
)
cell_cycle_genes <- unique(c(
  "Mki67", "Top2a", "Pcna", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6", "Mcm7",
  "Ccna2", "Ccnb1", "Ccnb2", "Ccnd1", "Ccne1", "Cdk1", "Cdk2", "Cdk4",
  "Birc5", "Ube2c", "Cenpf", "Nusap1", "Tyms", "Tk1", "Hmgb2"
))

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

scale_or_zero <- function(x, center, scale) {
  z <- (x - center) / scale
  z[!is.finite(z)] <- 0
  z
}

score_signature <- function(signature, deploy_train, deploy_heldout, train_meta) {
  if (nrow(signature) == 0) {
    return(list(raw = NA_real_, relative = NA_real_, coverage = 0, weighted_coverage = 0,
                young_center = NA_real_, old_center = NA_real_, denom = NA_real_))
  }
  genes <- intersect(signature$gene, rownames(deploy_train))
  genes <- intersect(genes, rownames(deploy_heldout))
  sig <- signature[match(genes, signature$gene), , drop = FALSE]
  if (nrow(sig) == 0) {
    return(list(raw = NA_real_, relative = NA_real_, coverage = 0, weighted_coverage = 0,
                young_center = NA_real_, old_center = NA_real_, denom = NA_real_))
  }
  mu <- rowMeans(deploy_train[genes, , drop = FALSE])
  sdv <- apply(deploy_train[genes, , drop = FALSE], 1, sd)
  usable <- is.finite(sdv) & sdv > 0 & is.finite(mu)
  genes <- genes[usable]
  sig <- sig[usable, , drop = FALSE]
  mu <- mu[usable]
  sdv <- sdv[usable]
  if (length(genes) == 0) {
    return(list(raw = NA_real_, relative = NA_real_, coverage = 0, weighted_coverage = 0,
                young_center = NA_real_, old_center = NA_real_, denom = NA_real_))
  }
  z_train <- sweep(sweep(deploy_train[genes, , drop = FALSE], 1, mu, "-"), 1, sdv, "/")
  z_held <- scale_or_zero(as.numeric(deploy_heldout[genes, 1]), mu, sdv)
  names(z_held) <- genes
  w <- sig$weight
  w[!is.finite(w) | w <= 0] <- min(w[is.finite(w) & w > 0], na.rm = TRUE)
  if (!any(is.finite(w) & w > 0)) w <- rep(1, nrow(sig))
  module_score <- function(zvec, module) {
    idx <- which(sig$module == module)
    if (length(idx) == 0) return(NA_real_)
    weighted.mean(zvec[idx], w[idx], na.rm = TRUE)
  }
  train_raw <- vapply(seq_len(ncol(z_train)), function(j) {
    module_score(z_train[, j], "young_high") - module_score(z_train[, j], "old_high")
  }, numeric(1))
  names(train_raw) <- colnames(z_train)
  held_raw <- module_score(z_held, "young_high") - module_score(z_held, "old_high")
  young_center <- median(train_raw[train_meta$age_group == "Young"], na.rm = TRUE)
  old_center <- median(train_raw[train_meta$age_group == "Old"], na.rm = TRUE)
  denom <- young_center - old_center
  relative <- if (is.finite(denom) && abs(denom) > epsilon) {
    (held_raw - old_center) / denom
  } else {
    NA_real_
  }
  total_weight <- sum(signature$weight, na.rm = TRUE)
  covered_weight <- sum(sig$weight, na.rm = TRUE)
  list(
    raw = held_raw,
    relative = relative,
    coverage = nrow(sig) / nrow(signature),
    weighted_coverage = covered_weight / total_weight,
    young_center = young_center,
    old_center = old_center,
    denom = denom
  )
}

build_signature <- function(ranked, per_direction) {
  young <- ranked[ranked$module == "young_high" & ranked$reliability_pass, , drop = FALSE]
  old <- ranked[ranked$module == "old_high" & ranked$reliability_pass, , drop = FALSE]
  young <- young[order(-young$rank_score, young$gene), , drop = FALSE]
  old <- old[order(-old$rank_score, old$gene), , drop = FALSE]
  sig <- rbind(head(young, per_direction), head(old, per_direction))
  sig$signature_target_per_direction <- per_direction
  sig
}

gene_flags <- function(genes) {
  data.frame(
    gene = genes,
    mitochondrial_flag = grepl("^mt-|^Mt-", genes),
    ribosomal_flag = grepl("^Rp[sl]", genes),
    cell_cycle_flag = genes %in% cell_cycle_genes,
    sex_linked_flag = genes %in% sex_linked_genes,
    stringsAsFactors = FALSE
  )
}

technical_burden <- function(sig, expression, train_meta) {
  if (nrow(sig) == 0) {
    return(data.frame(mean_abs_cor_library = NA_real_, mean_abs_cor_detected = NA_real_, mean_abs_cor_cell_count = NA_real_))
  }
  genes <- intersect(sig$gene, rownames(expression))
  if (length(genes) == 0) {
    return(data.frame(mean_abs_cor_library = NA_real_, mean_abs_cor_detected = NA_real_, mean_abs_cor_cell_count = NA_real_))
  }
  cors <- t(vapply(genes, function(g) {
    x <- as.numeric(expression[g, train_meta$mouse])
    c(
      abs(safe_cor(x, train_meta$pseudobulk_library_size)),
      abs(safe_cor(x, train_meta$pseudobulk_detected_genes)),
      abs(safe_cor(x, train_meta$n_cells))
    )
  }, numeric(3)))
  data.frame(
    mean_abs_cor_library = mean(cors[, 1], na.rm = TRUE),
    mean_abs_cor_detected = mean(cors[, 2], na.rm = TRUE),
    mean_abs_cor_cell_count = mean(cors[, 3], na.rm = TRUE)
  )
}

message("Reading FACS v2 pseudobulk inputs")
counts <- readRDS(counts_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(colnames(counts), metadata$mouse)) {
  stop("Count columns do not match metadata mouse order")
}
metadata$age_group <- factor(metadata$age_group, levels = c("Young", "Old"))
metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
metadata$age_sex_group <- interaction(metadata$sex, metadata$age_group, sep = "_", drop = TRUE)
metadata$young_label <- as.integer(metadata$age_group == "Young")

fold_score_rows <- list()
fold_diag_rows <- list()
fold_sig_rows <- list()
fold_gene_effect_rows <- list()
fold_pc_rows <- list()
fold_weak_rows <- list()

for (heldout in metadata$mouse) {
  message("Fold held out: ", heldout)
  held_idx <- which(metadata$mouse == heldout)
  train_idx <- setdiff(seq_len(nrow(metadata)), held_idx)
  train_meta <- metadata[train_idx, , drop = FALSE]
  held_meta <- metadata[held_idx, , drop = FALSE]
  train_meta$age_group <- factor(train_meta$age_group, levels = c("Young", "Old"))
  train_meta$sex <- factor(train_meta$sex, levels = c("female", "male"))
  train_meta$age_sex_group <- factor(
    interaction(train_meta$sex, train_meta$age_group, sep = "_", drop = TRUE),
    levels = c("female_Young", "male_Young", "female_Old", "male_Old")
  )
  group_counts <- table(train_meta$age_sex_group)
  expected_groups <- c("female_Young", "male_Young", "female_Old", "male_Old")
  support <- setNames(as.integer(group_counts[expected_groups]), paste0("train_n_", expected_groups))
  support[is.na(support)] <- 0L
  weak_support <- any(support < 2)

  train_counts_raw <- counts[, train_meta$mouse, drop = FALSE]
  held_counts_raw <- counts[, held_meta$mouse, drop = FALSE]
  train_lib <- colSums(train_counts_raw)
  train_cpm <- t(t(train_counts_raw) / train_lib) * 1e6
  keep <- rowSums(train_cpm > cpm_threshold) >= min_mice_expressed
  train_counts <- train_counts_raw[keep, , drop = FALSE]
  held_counts <- held_counts_raw[keep, , drop = FALSE]

  dge <- DGEList(counts = train_counts)
  dge <- calcNormFactors(dge, method = "TMM")
  design_age_only <- model.matrix(~ age_group, data = train_meta)
  design_factorial <- model.matrix(~ 0 + age_sex_group, data = train_meta)
  design_rank <- qr(design_factorial)$rank
  design_full_rank <- design_rank == ncol(design_factorial)
  residual_df <- nrow(design_factorial) - design_rank

  voom_age <- voom(dge, design_age_only, plot = FALSE)
  fit_age <- eBayes(lmFit(voom_age, design_age_only))
  age_tab <- topTable(fit_age, coef = "age_groupOld", number = Inf, sort.by = "none")

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

  genes <- rownames(train_counts)
  rownames(age_tab) <- rownames(female_tab) <- rownames(male_tab) <- rownames(common_tab) <- rownames(interaction_tab) <- genes
  ranked <- data.frame(
    gene = genes,
    age_only_logFC = age_tab[genes, "logFC"],
    age_only_t = age_tab[genes, "t"],
    age_only_adj.P.Val = age_tab[genes, "adj.P.Val"],
    female_logFC = female_tab[genes, "logFC"],
    female_t = female_tab[genes, "t"],
    male_logFC = male_tab[genes, "logFC"],
    male_t = male_tab[genes, "t"],
    common_logFC = common_tab[genes, "logFC"],
    common_t = common_tab[genes, "t"],
    common_adj.P.Val = common_tab[genes, "adj.P.Val"],
    interaction_logFC = interaction_tab[genes, "logFC"],
    interaction_t = interaction_tab[genes, "t"],
    interaction_adj.P.Val = interaction_tab[genes, "adj.P.Val"],
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
  ranked$interaction_penalty <- 1 / (1 + ranked$interaction_magnitude_ratio)
  ranked$reliability_pass <- with(
    ranked,
    sex_direction_concordant &
      age_only_common_concordant &
      !sex_linked_flag &
      is.finite(common_logFC) &
      is.finite(common_t)
  )
  ranked$module <- ifelse(ranked$common_logFC < 0, "young_high", ifelse(ranked$common_logFC > 0, "old_high", "neutral"))
  ranked$weight <- abs(ranked$common_logFC) * abs(ranked$common_t) * ranked$interaction_penalty
  ranked$rank_score <- ranked$weight * as.numeric(ranked$reliability_pass)

  deploy_train <- log2(t(t(train_counts) / colSums(train_counts)) * 1e6 + 1)
  deploy_held <- log2(t(t(held_counts) / colSums(held_counts)) * 1e6 + 1)
  stopifnot(identical(colnames(deploy_train), train_meta$mouse))
  training_expression_unchanged <- TRUE

  age_ranked <- ranked
  age_ranked$module <- ifelse(age_ranked$age_only_logFC < 0, "young_high", ifelse(age_ranked$age_only_logFC > 0, "old_high", "neutral"))
  age_ranked$reliability_pass <- with(age_ranked, !sex_linked_flag & age_only_common_concordant & is.finite(age_only_logFC) & is.finite(age_only_t))
  age_ranked$weight <- abs(age_ranked$age_only_logFC) * abs(age_ranked$age_only_t)
  age_ranked$rank_score <- age_ranked$weight * as.numeric(age_ranked$reliability_pass)

  train_logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  gene_var <- apply(train_logcpm, 1, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(top_variable_n, length(gene_var)))]
  pca <- prcomp(t(train_logcpm[top_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
  pc_scores <- data.frame(
    mouse = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  pc_scores <- merge(pc_scores, train_meta, by = "mouse", sort = FALSE)
  fold_pc_rows[[heldout]] <- data.frame(
    heldout_mouse = heldout,
    pca_gene_count = length(top_genes),
    pc1_cor_library = safe_cor(pc_scores$PC1, pc_scores$pseudobulk_library_size),
    pc1_cor_detected = safe_cor(pc_scores$PC1, pc_scores$pseudobulk_detected_genes),
    pc1_cor_cell_count = safe_cor(pc_scores$PC1, pc_scores$n_cells),
    pc2_cor_library = safe_cor(pc_scores$PC2, pc_scores$pseudobulk_library_size),
    pc2_cor_detected = safe_cor(pc_scores$PC2, pc_scores$pseudobulk_detected_genes),
    pc2_cor_cell_count = safe_cor(pc_scores$PC2, pc_scores$n_cells),
    stringsAsFactors = FALSE
  )

  fold_gene_effect_rows[[heldout]] <- data.frame(
    heldout_mouse = heldout,
    gene = ranked$gene,
    eligible = ranked$reliability_pass,
    module = ranked$module,
    common_direction = ranked$common_direction,
    female_direction = ranked$female_direction,
    male_direction = ranked$male_direction,
    age_only_direction = ranked$age_only_direction,
    common_logFC = ranked$common_logFC,
    female_logFC = ranked$female_logFC,
    male_logFC = ranked$male_logFC,
    interaction_penalty = ranked$interaction_penalty,
    low_depth_direction_checked = FALSE,
    stringsAsFactors = FALSE
  )

  fold_diag <- data.frame(
    heldout_mouse = heldout,
    heldout_age = held_meta$age,
    heldout_age_months = held_meta$age_months,
    heldout_age_group = as.character(held_meta$age_group),
    heldout_sex = as.character(held_meta$sex),
    fold_type = ifelse(weak_support, "weak_support", "standard_support"),
    train_n_mice = nrow(train_meta),
    train_n_young = sum(train_meta$age_group == "Young"),
    train_n_old = sum(train_meta$age_group == "Old"),
    t(support),
    design_full_rank = design_full_rank,
    residual_df = residual_df,
    weak_support_fold = weak_support,
    retained_genes = sum(keep),
    reliability_pass_genes = sum(ranked$reliability_pass),
    young_high_candidates = sum(ranked$reliability_pass & ranked$module == "young_high"),
    old_high_candidates = sum(ranked$reliability_pass & ranked$module == "old_high"),
    heldout_raw_library_size = held_meta$pseudobulk_library_size,
    heldout_detected_genes = held_meta$pseudobulk_detected_genes,
    heldout_cell_count = held_meta$n_cells,
    heldout_normalization_external = TRUE,
    training_expression_unchanged = training_expression_unchanged,
    stringsAsFactors = FALSE
  )

  for (version in names(signature_sizes)) {
    sig <- build_signature(ranked, signature_sizes[[version]])
    age_sig <- build_signature(age_ranked, signature_sizes[[version]])
    sig_flags <- gene_flags(sig$gene)
    sig$mitochondrial_flag <- sig_flags$mitochondrial_flag
    sig$ribosomal_flag <- sig_flags$ribosomal_flag
    sig$cell_cycle_flag <- sig_flags$cell_cycle_flag
    sig$sex_linked_flag <- sig_flags$sex_linked_flag
    sig$heldout_mouse <- heldout
    sig$signature_version <- version
    fold_sig_rows[[paste(heldout, version, sep = "_")]] <- sig

    sc <- score_signature(sig, deploy_train, deploy_held, train_meta)
    age_sc <- score_signature(age_sig, deploy_train, deploy_held, train_meta)
    tech <- technical_burden(sig, deploy_train, train_meta)
    fold_score_rows[[paste(heldout, version, sep = "_")]] <- cbind(
      data.frame(
        mouse = heldout,
        age = held_meta$age,
        age_months = held_meta$age_months,
        age_group = as.character(held_meta$age_group),
        sex = as.character(held_meta$sex),
        fold_type = ifelse(weak_support, "weak_support", "standard_support"),
        signature_version = version,
        signature_size = nrow(sig),
        young_high_n = sum(sig$module == "young_high"),
        old_high_n = sum(sig$module == "old_high"),
        predicted_raw_score = sc$raw,
        predicted_calibrated_score = sc$relative,
        age_only_predicted_raw_score = age_sc$raw,
        age_only_predicted_calibrated_score = age_sc$relative,
        gene_coverage = sc$coverage,
        weighted_gene_coverage = sc$weighted_coverage,
        young_reference_center = sc$young_center,
        old_reference_center = sc$old_center,
        calibration_denominator = sc$denom,
        mitochondrial_fraction = mean(sig$mitochondrial_flag),
        ribosomal_fraction = mean(sig$ribosomal_flag),
        cell_cycle_fraction = mean(sig$cell_cycle_flag),
        sex_linked_fraction = mean(sig$sex_linked_flag),
        heldout_raw_library_size = held_meta$pseudobulk_library_size,
        heldout_detected_genes = held_meta$pseudobulk_detected_genes,
        heldout_cell_count = held_meta$n_cells,
        stringsAsFactors = FALSE
      ),
      tech
    )
  }

  weak_stats <- aggregate(
    cbind(common_logFC, female_logFC, male_logFC, interaction_penalty) ~ 1,
    data = ranked[ranked$reliability_pass, ],
    FUN = function(x) c(median = median(x, na.rm = TRUE), q25 = quantile(x, 0.25, na.rm = TRUE), q75 = quantile(x, 0.75, na.rm = TRUE))
  )
  fold_weak_rows[[heldout]] <- data.frame(
    heldout_mouse = heldout,
    weak_support_fold = weak_support,
    reliability_pass_genes = sum(ranked$reliability_pass),
    median_common_logFC = median(ranked$common_logFC[ranked$reliability_pass], na.rm = TRUE),
    median_female_logFC = median(ranked$female_logFC[ranked$reliability_pass], na.rm = TRUE),
    median_male_logFC = median(ranked$male_logFC[ranked$reliability_pass], na.rm = TRUE),
    median_interaction_penalty = median(ranked$interaction_penalty[ranked$reliability_pass], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  fold_diag_rows[[heldout]] <- fold_diag
}

scores <- do.call(rbind, fold_score_rows)
diagnostics <- do.call(rbind, fold_diag_rows)
signatures <- do.call(rbind, fold_sig_rows)
gene_effects <- do.call(rbind, fold_gene_effect_rows)
pc_audit <- do.call(rbind, fold_pc_rows)
weak_audit <- do.call(rbind, fold_weak_rows)

write.csv(scores, file.path(validation_dir, "step11_nested_lomo_scores.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(validation_dir, "step11_nested_lomo_fold_diagnostics.csv"), row.names = FALSE)
write.csv(signatures, file.path(validation_dir, "step11_nested_lomo_fold_signatures.csv"), row.names = FALSE)
write.csv(gene_effects, file.path(validation_dir, "step11_nested_lomo_gene_effects_by_fold.csv"), row.names = FALSE)
write.csv(pc_audit, file.path(validation_dir, "step11_nested_lomo_fold_pc_technical_audit.csv"), row.names = FALSE)
write.csv(weak_audit, file.path(validation_dir, "step11_nested_lomo_weak_support_audit.csv"), row.names = FALSE)

message("Summarizing nested LOMO feasibility")
summary_rows <- list()
for (version in names(signature_sizes)) {
  sv <- scores[scores$signature_version == version, , drop = FALSE]
  summary_rows[[version]] <- data.frame(
    signature_version = version,
    n_folds = nrow(sv),
    zero_signature_folds = sum(sv$signature_size == 0),
    median_signature_size = median(sv$signature_size),
    all_spearman_age = safe_cor(sv$predicted_calibrated_score, sv$age_months),
    old_only_spearman_age = safe_cor(sv$predicted_calibrated_score[sv$age_group == "Old"], sv$age_months[sv$age_group == "Old"]),
    young_old_auc = auc_rank(as.integer(sv$age_group == "Young"), sv$predicted_calibrated_score),
    young_median = median(sv$predicted_calibrated_score[sv$age_group == "Young"], na.rm = TRUE),
    old_median = median(sv$predicted_calibrated_score[sv$age_group == "Old"], na.rm = TRUE),
    score_cor_library = safe_cor(sv$predicted_calibrated_score, sv$heldout_raw_library_size),
    score_cor_detected = safe_cor(sv$predicted_calibrated_score, sv$heldout_detected_genes),
    score_cor_cell_count = safe_cor(sv$predicted_calibrated_score, sv$heldout_cell_count),
    median_gene_coverage = median(sv$gene_coverage, na.rm = TRUE),
    median_weighted_gene_coverage = median(sv$weighted_gene_coverage, na.rm = TRUE),
    median_ribosomal_fraction = median(sv$ribosomal_fraction, na.rm = TRUE),
    median_cell_cycle_fraction = median(sv$cell_cycle_fraction, na.rm = TRUE),
    median_gene_library_burden = median(sv$mean_abs_cor_library, na.rm = TRUE),
    median_gene_detected_burden = median(sv$mean_abs_cor_detected, na.rm = TRUE),
    median_gene_cell_count_burden = median(sv$mean_abs_cor_cell_count, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(validation_dir, "step11_nested_lomo_feasibility_summary.csv"), row.names = FALSE)

signature_sets <- split(signatures$gene, paste(signatures$signature_version, signatures$heldout_mouse, sep = "::"))
jaccard_rows <- list()
for (version in names(signature_sizes)) {
  keys <- grep(paste0("^", version, "::"), names(signature_sets), value = TRUE)
  pairs <- combn(keys, 2, simplify = FALSE)
  jaccard_rows[[version]] <- do.call(rbind, lapply(pairs, function(pair) {
    a <- unique(signature_sets[[pair[1]]])
    b <- unique(signature_sets[[pair[2]]])
    data.frame(
      signature_version = version,
      fold_a = sub(".*::", "", pair[1]),
      fold_b = sub(".*::", "", pair[2]),
      jaccard = length(intersect(a, b)) / length(union(a, b)),
      stringsAsFactors = FALSE
    )
  }))
}
jaccard <- do.call(rbind, jaccard_rows)
write.csv(jaccard, file.path(validation_dir, "step11_nested_lomo_signature_jaccard.csv"), row.names = FALSE)

selection_freq <- aggregate(
  heldout_mouse ~ signature_version + gene + module,
  data = signatures,
  FUN = length
)
colnames(selection_freq)[colnames(selection_freq) == "heldout_mouse"] <- "selected_folds"
selection_freq$selection_frequency <- selection_freq$selected_folds / length(unique(metadata$mouse))
write.csv(selection_freq[order(selection_freq$signature_version, -selection_freq$selection_frequency, selection_freq$gene), ],
          file.path(validation_dir, "step11_nested_lomo_gene_selection_frequency.csv"), row.names = FALSE)

png(file.path(validation_dir, "step11_nested_lomo_scores_by_age.png"), width = 1200, height = 850, res = 150)
plot(
  NA,
  xlim = range(scores$age_months),
  ylim = range(scores$predicted_calibrated_score, na.rm = TRUE),
  xlab = "Age months",
  ylab = "Nested LOMO calibrated score",
  main = "FACS v2 nested LOMO feasibility scores"
)
cols <- c(Small = "#7A5195", Medium = "#2878B5", Large = "#D95F02")
pchs <- ifelse(scores$fold_type == "weak_support", 17, 16)
for (version in names(signature_sizes)) {
  sv <- scores[scores$signature_version == version, ]
  points(sv$age_months + match(version, names(signature_sizes)) * 0.22 - 0.44, sv$predicted_calibrated_score,
         col = cols[[version]], pch = pchs[scores$signature_version == version])
}
legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n")
dev.off()

lowest_depth_mouse <- metadata$mouse[which.min(metadata$pseudobulk_library_size)]
removal_summary <- data.frame()
for (version in names(signature_sizes)) {
  sv <- scores[scores$signature_version == version, ]
  sv2 <- sv[sv$mouse != lowest_depth_mouse, ]
  removal_summary <- rbind(removal_summary, data.frame(
    signature_version = version,
    removed_mouse = lowest_depth_mouse,
    all_spearman_age_after_removal = safe_cor(sv2$predicted_calibrated_score, sv2$age_months),
    old_only_spearman_age_after_removal = safe_cor(sv2$predicted_calibrated_score[sv2$age_group == "Old"], sv2$age_months[sv2$age_group == "Old"]),
    auc_after_removal = auc_rank(as.integer(sv2$age_group == "Young"), sv2$predicted_calibrated_score),
    score_cor_library_after_removal = safe_cor(sv2$predicted_calibrated_score, sv2$heldout_raw_library_size),
    stringsAsFactors = FALSE
  ))
}
write.csv(removal_summary, file.path(validation_dir, "step11_lowest_depth_removal_sensitivity.csv"), row.names = FALSE)

report_lines <- c(
  "# FACS v2 Step 11: Fully Nested LOMO Feasibility Audit",
  "",
  "## Scope",
  "",
  "This is the first fully nested LOMO run for feasibility auditing. It is not the final model-selection or permutation stage.",
  "",
  "Every outer fold repeats training-side gene filtering, TMM, voom, age-only DE, factorial DE, contrasts, sex-linked exclusion, reliability audit, ranking, signature selection, weights, frozen deployable normalization, and held-out scoring.",
  "",
  "No files under `data_facs` were modified.",
  "",
  "## Main Feasibility Summary",
  "",
  paste(capture.output(print(summary_table, row.names = FALSE)), collapse = "\n"),
  "",
  "## Fold Diagnostics",
  "",
  paste0("- Folds: ", nrow(diagnostics)),
  paste0("- Full-rank factorial folds: ", sum(diagnostics$design_full_rank), " / ", nrow(diagnostics)),
  paste0("- Weak-support folds: ", sum(diagnostics$weak_support_fold), " / ", nrow(diagnostics)),
  paste0("- Zero-signature folds across all versions: ", sum(scores$signature_size == 0)),
  paste0("- Lowest-depth mouse for sensitivity: ", lowest_depth_mouse),
  "",
  "## Technical Axis Audit",
  "",
  "Per-fold PC1/PC2 technical correlations are recorded in `step11_nested_lomo_fold_pc_technical_audit.csv`.",
  "Per-score held-out technical associations are summarized above and retained per mouse in `step11_nested_lomo_scores.csv`.",
  "",
  "## Weak-Support Fold Handling",
  "",
  "Weak-support folds are not removed. They are flagged for separate interpretation because one age-by-sex training cell has only one mouse.",
  "",
  "## Stability Gate Candidates",
  "",
  "This run records gene-level fold effects and selection frequency. Final primary candidates should require high common-direction consistency, male direction stability, no systematic female reversal, age-only/common agreement, and adequate fold-selection frequency.",
  "",
  "## Outputs",
  "",
  "- `outputs/facs_v2/validation/step11_nested_lomo_scores.csv`",
  "- `outputs/facs_v2/validation/step11_nested_lomo_fold_diagnostics.csv`",
  "- `outputs/facs_v2/validation/step11_nested_lomo_fold_signatures.csv`",
  "- `outputs/facs_v2/validation/step11_nested_lomo_gene_effects_by_fold.csv`",
  "- `outputs/facs_v2/validation/step11_nested_lomo_fold_pc_technical_audit.csv`",
  "- `outputs/facs_v2/validation/step11_nested_lomo_signature_jaccard.csv`",
  "- `outputs/facs_v2/validation/step11_nested_lomo_gene_selection_frequency.csv`",
  "- `outputs/facs_v2/validation/step11_lowest_depth_removal_sensitivity.csv`"
)
writeLines(report_lines, file.path(validation_dir, "step11_nested_lomo_feasibility_report.md"))
writeLines(capture.output(sessionInfo()), file.path(validation_dir, "sessionInfo_step11.txt"))

message("Done")
message(sprintf("Report: %s", file.path(validation_dir, "step11_nested_lomo_feasibility_report.md")))
