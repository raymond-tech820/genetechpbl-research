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
model_comp_dir <- file.path(out_root, "model_comparison")
out_dir <- file.path(out_root, "structural_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_all <- readRDS(file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.rds"))
metadata_all <- read.csv(file.path(processed_dir, "facs_v2_limb_msc_mouse_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(colnames(counts_all), metadata_all$mouse))

cpm_threshold <- 1
min_mice_expressed <- 2
epsilon <- 1e-6
medium_per_direction <- 50
sex_linked_genes <- c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Rps4y1", "Rps4y2", "Sry", "Jarid1d")

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

prep_meta <- function(meta) {
  meta$age_group <- ifelse(meta$age == "3m", "Young", "Old")
  meta$age_group <- factor(meta$age_group, levels = c("Young", "Old"))
  meta$sex <- factor(meta$sex, levels = c("female", "male"))
  meta$age_sex_group <- factor(interaction(meta$sex, meta$age_group, sep = "_", drop = TRUE),
                               levels = c("female_Young", "male_Young", "female_Old", "male_Old"))
  meta
}

rank_factorial <- function(train_counts, train_meta) {
  dge <- DGEList(counts = train_counts)
  dge <- calcNormFactors(dge, method = "TMM")
  design_age <- model.matrix(~ age_group, data = train_meta)
  design_fact <- model.matrix(~ 0 + age_sex_group, data = train_meta)
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
  female <- topTable(fit_con, coef = "female_age", number = Inf, sort.by = "none")
  male <- topTable(fit_con, coef = "male_age", number = Inf, sort.by = "none")
  common <- topTable(fit_con, coef = "common_age", number = Inf, sort.by = "none")
  inter <- topTable(fit_con, coef = "interaction_age_by_sex", number = Inf, sort.by = "none")
  genes <- rownames(train_counts)
  rownames(age_tab) <- rownames(female) <- rownames(male) <- rownames(common) <- rownames(inter) <- genes
  r <- data.frame(
    gene = genes,
    age_only_logFC = age_tab[genes, "logFC"],
    age_only_t = age_tab[genes, "t"],
    female_logFC = female[genes, "logFC"],
    female_t = female[genes, "t"],
    male_logFC = male[genes, "logFC"],
    male_t = male[genes, "t"],
    common_logFC = common[genes, "logFC"],
    common_t = common[genes, "t"],
    interaction_logFC = inter[genes, "logFC"],
    stringsAsFactors = FALSE
  )
  r$sex_linked_flag <- r$gene %in% sex_linked_genes
  r$female_direction <- sign(r$female_logFC)
  r$male_direction <- sign(r$male_logFC)
  r$common_direction <- sign(r$common_logFC)
  r$age_only_direction <- sign(r$age_only_logFC)
  r$sex_direction_concordant <- r$female_direction == r$male_direction & r$female_direction != 0
  r$age_only_common_concordant <- r$age_only_direction == r$common_direction & r$common_direction != 0
  r$interaction_ratio <- abs(r$interaction_logFC) / (abs(r$female_logFC) + abs(r$male_logFC) + epsilon)
  r$interaction_penalty <- 1 / (1 + r$interaction_ratio)
  r$reliability_pass <- with(r, sex_direction_concordant & age_only_common_concordant & !sex_linked_flag & is.finite(common_logFC) & is.finite(common_t))
  r$module <- ifelse(r$common_logFC < 0, "young_high", ifelse(r$common_logFC > 0, "old_high", "neutral"))
  r$weight <- abs(r$common_logFC) * abs(r$common_t) * r$interaction_penalty
  r$rank_score <- r$weight * as.numeric(r$reliability_pass)
  r
}

build_medium <- function(ranked) {
  r <- ranked[ranked$reliability_pass, , drop = FALSE]
  y <- r[r$module == "young_high", , drop = FALSE]
  o <- r[r$module == "old_high", , drop = FALSE]
  y <- y[order(-y$rank_score, y$gene), , drop = FALSE]
  o <- o[order(-o$rank_score, o$gene), , drop = FALSE]
  rbind(head(y, medium_per_direction), head(o, medium_per_direction))
}

filter_counts <- function(counts, train_mice) {
  raw <- counts[, train_mice, drop = FALSE]
  keep <- rowSums(t(t(raw) / colSums(raw)) * 1e6 > cpm_threshold) >= min_mice_expressed
  keep
}

score_sig <- function(sig, train_counts, held_counts, train_meta) {
  train_expr <- log2(t(t(train_counts) / colSums(train_counts)) * 1e6 + 1)
  held_expr <- log2(t(t(held_counts) / colSums(held_counts)) * 1e6 + 1)
  genes <- intersect(sig$gene, rownames(train_expr))
  genes <- intersect(genes, rownames(held_expr))
  sig <- sig[match(genes, sig$gene), , drop = FALSE]
  if (length(genes) == 0 || nrow(sig) == 0) return(c(raw = NA, relative = NA, denom = NA, coverage = 0))
  mu <- rowMeans(train_expr[genes, , drop = FALSE])
  sdv <- apply(train_expr[genes, , drop = FALSE], 1, sd)
  usable <- is.finite(mu) & is.finite(sdv) & sdv > 0
  genes <- genes[usable]; sig <- sig[usable, , drop = FALSE]; mu <- mu[usable]; sdv <- sdv[usable]
  ztrain <- sweep(sweep(train_expr[genes, , drop = FALSE], 1, mu, "-"), 1, sdv, "/")
  zheld <- (as.numeric(held_expr[genes, 1]) - mu) / sdv
  w <- sig$weight
  w[!is.finite(w) | w <= 0] <- 1
  ms <- function(z, module) {
    idx <- which(sig$module == module)
    if (!length(idx)) return(NA_real_)
    weighted.mean(z[idx], w[idx], na.rm = TRUE)
  }
  train_raw <- vapply(seq_len(ncol(ztrain)), function(i) ms(ztrain[, i], "young_high") - ms(ztrain[, i], "old_high"), numeric(1))
  held_raw <- ms(zheld, "young_high") - ms(zheld, "old_high")
  yc <- median(train_raw[train_meta$age_group == "Young"], na.rm = TRUE)
  oc <- median(train_raw[train_meta$age_group == "Old"], na.rm = TRUE)
  denom <- yc - oc
  rel <- if (is.finite(denom) && abs(denom) > epsilon) (held_raw - oc) / denom else NA_real_
  c(raw = held_raw, relative = rel, denom = denom, coverage = nrow(sig) / 100)
}

run_3m18_lomo <- function(counts, meta) {
  rows <- list(); sigs <- list(); effects <- list(); diags <- list()
  for (heldout in meta$mouse) {
    message("3m-vs-18m held out: ", heldout)
    train_meta <- meta[meta$mouse != heldout, , drop = FALSE]
    held_meta <- meta[meta$mouse == heldout, , drop = FALSE]
    support <- table(train_meta$age_sex_group)
    full_rank <- qr(model.matrix(~ 0 + age_sex_group, data = train_meta))$rank == 4
    keep <- filter_counts(counts, train_meta$mouse)
    train_counts <- counts[keep, train_meta$mouse, drop = FALSE]
    held_counts <- counts[keep, heldout, drop = FALSE]
    ranked <- rank_factorial(train_counts, train_meta)
    sig <- build_medium(ranked)
    sc <- score_sig(sig, train_counts, held_counts, train_meta)
    rows[[heldout]] <- data.frame(
      mouse = heldout,
      age = held_meta$age,
      age_months = held_meta$age_months,
      age_group = as.character(held_meta$age_group),
      sex = as.character(held_meta$sex),
      score = sc[["relative"]],
      raw_score = sc[["raw"]],
      signature_size = nrow(sig),
      young_high_n = sum(sig$module == "young_high"),
      old_high_n = sum(sig$module == "old_high"),
      library_size = held_meta$pseudobulk_library_size,
      detected_genes = held_meta$pseudobulk_detected_genes,
      cell_count = held_meta$n_cells,
      calibration_denominator = sc[["denom"]],
      gene_coverage = sc[["coverage"]],
      design_full_rank = full_rank,
      weak_support_fold = any(support < 2),
      stringsAsFactors = FALSE
    )
    sig$heldout_mouse <- heldout
    sigs[[heldout]] <- sig
    effects[[heldout]] <- transform(ranked, heldout_mouse = heldout)
    diags[[heldout]] <- data.frame(heldout_mouse = heldout, retained_genes = sum(keep), reliability_pass = sum(ranked$reliability_pass), design_full_rank = full_rank)
  }
  list(scores = do.call(rbind, rows), signatures = do.call(rbind, sigs), effects = do.call(rbind, effects), diagnostics = do.call(rbind, diags))
}

summarize_scores <- function(scores) {
  data.frame(
    n_mice = nrow(scores),
    auc = auc_rank(as.integer(scores$age_group == "Young"), scores$score),
    all_rho = safe_cor(scores$score, scores$age_months),
    library_rho = safe_cor(scores$score, scores$library_size),
    detected_rho = safe_cor(scores$score, scores$detected_genes),
    cell_count_rho = safe_cor(scores$score, scores$cell_count),
    young_minus_old_median = median(scores$score[scores$age_group == "Young"], na.rm = TRUE) - median(scores$score[scores$age_group == "Old"], na.rm = TRUE),
    zero_signature_folds = sum(scores$signature_size == 0),
    full_rank_folds = sum(scores$design_full_rank),
    weak_support_folds = sum(scores$weak_support_fold)
  )
}

jaccard_median <- function(sigs) {
  sets <- split(sigs$gene, sigs$heldout_mouse)
  pairs <- combn(names(sets), 2, simplify = FALSE)
  median(vapply(pairs, function(p) length(intersect(sets[[p[1]]], sets[[p[2]]])) / length(union(sets[[p[1]]], sets[[p[2]]])), numeric(1)))
}

message("Preparing 3m and 18m subset")
meta_3_18 <- metadata_all[metadata_all$age %in% c("3m", "18m"), , drop = FALSE]
meta_3_18 <- prep_meta(meta_3_18)
counts_3_18 <- counts_all[, meta_3_18$mouse, drop = FALSE]
res <- run_3m18_lomo(counts_3_18, meta_3_18)
summary_3_18 <- summarize_scores(res$scores)
summary_3_18$median_jaccard <- jaccard_median(res$signatures)

message("Training frozen 3m/18m model and applying to 24m")
keep_full <- filter_counts(counts_all, meta_3_18$mouse)
train_counts <- counts_all[keep_full, meta_3_18$mouse, drop = FALSE]
ranked_full <- rank_factorial(train_counts, meta_3_18)
sig_full <- build_medium(ranked_full)
apply_mice <- metadata_all$mouse[metadata_all$age == "24m"]
extension_rows <- list()
for (m in c(meta_3_18$mouse, apply_mice)) {
  held_counts <- counts_all[keep_full, m, drop = FALSE]
  held_meta <- metadata_all[metadata_all$mouse == m, , drop = FALSE]
  sc <- score_sig(sig_full, train_counts, held_counts, meta_3_18)
  extension_rows[[m]] <- data.frame(
    mouse = m,
    age = held_meta$age,
    age_months = held_meta$age_months,
    age_group_training_definition = ifelse(held_meta$age == "3m", "Young", ifelse(held_meta$age == "18m", "Old_train", "Extension_24m")),
    sex = held_meta$sex,
    score = sc[["relative"]],
    raw_score = sc[["raw"]],
    library_size = held_meta$pseudobulk_library_size,
    detected_genes = held_meta$pseudobulk_detected_genes,
    cell_count = held_meta$n_cells,
    stringsAsFactors = FALSE
  )
}
extension_scores <- do.call(rbind, extension_rows)
extension_summary <- aggregate(score ~ age + age_months + sex, extension_scores, median)
male_ext <- extension_scores[extension_scores$sex == "male", ]
male_ext_summary <- data.frame(
  male_age_spearman = safe_cor(male_ext$score, male_ext$age_months),
  median_3m_male = median(male_ext$score[male_ext$age == "3m"]),
  median_18m_male = median(male_ext$score[male_ext$age == "18m"]),
  median_24m_male = median(male_ext$score[male_ext$age == "24m"]),
  expected_3_gt_18_gt_24 = median(male_ext$score[male_ext$age == "3m"]) > median(male_ext$score[male_ext$age == "18m"]) &&
    median(male_ext$score[male_ext$age == "18m"]) > median(male_ext$score[male_ext$age == "24m"])
)

pooled_effect_path <- file.path(out_root, "de", "step07_10_full_data_factorial_de_and_ranking.csv")
pooled_effect <- if (file.exists(pooled_effect_path)) read.csv(pooled_effect_path, stringsAsFactors = FALSE) else data.frame()
effect_compare <- data.frame()
if (nrow(pooled_effect)) {
  full_3_18_effect <- ranked_full[, c("gene", "common_logFC", "female_logFC", "male_logFC")]
  names(full_3_18_effect)[-1] <- paste0("three18_", names(full_3_18_effect)[-1])
  pooled_small <- pooled_effect[, c("gene", "common_logFC", "female_logFC", "male_logFC")]
  names(pooled_small)[-1] <- paste0("pooled_", names(pooled_small)[-1])
  merged <- merge(full_3_18_effect, pooled_small, by = "gene")
  effect_compare <- data.frame(
    common_effect_spearman = safe_cor(merged$three18_common_logFC, merged$pooled_common_logFC),
    female_effect_spearman = safe_cor(merged$three18_female_logFC, merged$pooled_female_logFC),
    male_effect_spearman = safe_cor(merged$three18_male_logFC, merged$pooled_male_logFC),
    compared_genes = nrow(merged)
  )
}

model_comp_sig_path <- file.path(model_comp_dir, "model_comparison_fold_signatures.csv")
signature_overlap <- data.frame()
if (file.exists(model_comp_sig_path)) {
  pooled_sigs <- read.csv(model_comp_sig_path, stringsAsFactors = FALSE)
  pooled_medium <- pooled_sigs[pooled_sigs$model == "factorial_medium_original", ]
  overlap_by_fold <- merge(
    aggregate(gene ~ heldout_mouse, res$signatures, function(x) paste(unique(x), collapse = ";")),
    aggregate(gene ~ heldout_mouse, pooled_medium, function(x) paste(unique(x), collapse = ";")),
    by = "heldout_mouse",
    suffixes = c("_three18", "_pooled")
  )
  signature_overlap <- do.call(rbind, lapply(seq_len(nrow(overlap_by_fold)), function(i) {
    a <- strsplit(overlap_by_fold$gene_three18[i], ";", fixed = TRUE)[[1]]
    b <- strsplit(overlap_by_fold$gene_pooled[i], ";", fixed = TRUE)[[1]]
    data.frame(heldout_mouse = overlap_by_fold$heldout_mouse[i], medium_signature_jaccard_vs_pooled = length(intersect(a, b)) / length(union(a, b)))
  }))
}

write.csv(res$scores, file.path(out_dir, "three_m_vs_18m_nested_lomo_scores.csv"), row.names = FALSE)
write.csv(res$signatures, file.path(out_dir, "three_m_vs_18m_fold_signatures.csv"), row.names = FALSE)
write.csv(res$effects, file.path(out_dir, "three_m_vs_18m_gene_effects_by_fold.csv"), row.names = FALSE)
write.csv(res$diagnostics, file.path(out_dir, "three_m_vs_18m_fold_diagnostics.csv"), row.names = FALSE)
write.csv(summary_3_18, file.path(out_dir, "three_m_vs_18m_summary.csv"), row.names = FALSE)
write.csv(extension_scores, file.path(out_dir, "three18_frozen_model_24m_extension_scores.csv"), row.names = FALSE)
write.csv(extension_summary, file.path(out_dir, "three18_frozen_model_24m_extension_age_summary.csv"), row.names = FALSE)
write.csv(male_ext_summary, file.path(out_dir, "three18_frozen_model_24m_extension_male_summary.csv"), row.names = FALSE)
write.csv(effect_compare, file.path(out_dir, "three_m_vs_18m_effect_correlation_vs_pooled.csv"), row.names = FALSE)
write.csv(signature_overlap, file.path(out_dir, "three_m_vs_18m_signature_overlap_vs_pooled.csv"), row.names = FALSE)

report <- c(
  "# FACS v2 Step 15: 3m-vs-18m Structural Sensitivity and 24m Extension",
  "",
  "## Scope",
  "",
  "This step removes the male-only 24m group from training to test whether pooled-old design drives instability. The frozen 3m/18m model is then applied to 24m male mice without refitting or recalibration.",
  "",
  "No files under `data_facs` were modified.",
  "",
  "## 3m-vs-18m Nested LOMO Summary",
  "",
  paste(capture.output(print(summary_3_18, row.names = FALSE)), collapse = "\n"),
  "",
  "## Effect Correlation vs Pooled 18/24m Model",
  "",
  paste(capture.output(print(effect_compare, row.names = FALSE)), collapse = "\n"),
  "",
  "## Frozen 3m/18m -> 24m Male Extension",
  "",
  paste(capture.output(print(extension_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Male Age Trajectory Check",
  "",
  paste(capture.output(print(male_ext_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Outputs",
  "",
  "- `three_m_vs_18m_nested_lomo_scores.csv`",
  "- `three_m_vs_18m_fold_signatures.csv`",
  "- `three_m_vs_18m_gene_effects_by_fold.csv`",
  "- `three_m_vs_18m_summary.csv`",
  "- `three18_frozen_model_24m_extension_scores.csv`",
  "- `three18_frozen_model_24m_extension_male_summary.csv`",
  "- `three_m_vs_18m_effect_correlation_vs_pooled.csv`",
  "- `three_m_vs_18m_signature_overlap_vs_pooled.csv`"
)
writeLines(report, file.path(out_dir, "three_m_vs_18m_sensitivity_and_24m_extension_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step15.txt"))

message("Done")
message(sprintf("Report: %s", file.path(out_dir, "three_m_vs_18m_sensitivity_and_24m_extension_report.md")))
