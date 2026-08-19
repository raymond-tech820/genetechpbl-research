#!/usr/bin/env Rscript

# FACS Youth Score v2 Step 18: train and export full-data frozen models.
# Primary: Factorial Medium Original.
# Comparators: Large Original, Stability-Selected, Medium Equal-Weight, Age-Only.
# Raw data under data_facs is not read or modified; inputs are processed pseudobulk outputs.

options(stringsAsFactors = FALSE)

root <- "/Users/strangenoah/Desktop/AI+X/Genetech/Youth_score"
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

set.seed(20260721)

out_root <- file.path(root, "outputs", "facs_v2")
processed_dir <- file.path(out_root, "processed")
model_comparison_dir <- file.path(out_root, "model_comparison")
bootstrap_dir <- file.path(out_root, "bootstrap")
permutation_dir <- file.path(out_root, "permutation", "formal_medium_999")
out_dir <- file.path(out_root, "final_models")
model_dir <- file.path(out_dir, "models")
parser_dir <- file.path(out_dir, "parser")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(parser_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.rds"))
metadata <- read.csv(file.path(processed_dir, "facs_v2_limb_msc_mouse_metadata.csv"), check.names = FALSE)
stopifnot(identical(colnames(counts), metadata$mouse))
metadata$age_group <- factor(metadata$age_group, levels = c("Young", "Old"))
metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
metadata$age_sex_group <- factor(interaction(metadata$sex, metadata$age_group, sep = "_", drop = TRUE), levels = c("female_Young", "male_Young", "female_Old", "male_Old"))

cpm_threshold <- 1
min_mice_expressed <- 2
epsilon <- 1e-6
sex_linked_genes <- c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Rps4y1", "Rps4y2", "Sry", "Jarid1d")

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

rbind_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  cols <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(df) {
    missing <- setdiff(cols, names(df))
    for (m in missing) df[[m]] <- NA
    df[, cols, drop = FALSE]
  })
  do.call(rbind, x)
}


prepare_full <- function(count_mat, meta) {
  raw <- count_mat[, meta$mouse, drop = FALSE]
  keep <- rowSums(t(t(raw) / colSums(raw)) * 1e6 > cpm_threshold) >= min_mice_expressed
  train_counts <- raw[keep, , drop = FALSE]
  dge <- DGEList(counts = train_counts)
  dge <- normLibSizes(dge, method = "TMM")
  list(counts = train_counts, dge = dge, keep = keep)
}

rank_factorial <- function(train_counts, train_meta) {
  train_meta$age_group <- factor(train_meta$age_group, levels = c("Young", "Old"))
  train_meta$sex <- factor(train_meta$sex, levels = c("female", "male"))
  train_meta$age_sex_group <- factor(interaction(train_meta$sex, train_meta$age_group, sep = "_", drop = TRUE), levels = c("female_Young", "male_Young", "female_Old", "male_Old"))
  design_age <- model.matrix(~ age_group, data = train_meta)
  design_fact <- model.matrix(~ 0 + age_sex_group, data = train_meta)
  if (qr(design_age)$rank < ncol(design_age) || qr(design_fact)$rank < ncol(design_fact)) stop("Full-data design is not full rank")

  dge <- DGEList(counts = train_counts)
  dge <- normLibSizes(dge, method = "TMM")
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
    age_only_p_value = age_tab[genes, "P.Value"],
    age_only_adj_p_value = age_tab[genes, "adj.P.Val"],
    female_logFC = female[genes, "logFC"],
    female_t = female[genes, "t"],
    male_logFC = male[genes, "logFC"],
    male_t = male[genes, "t"],
    common_logFC = common[genes, "logFC"],
    common_t = common[genes, "t"],
    common_p_value = common[genes, "P.Value"],
    common_adj_p_value = common[genes, "adj.P.Val"],
    interaction_logFC = inter[genes, "logFC"],
    interaction_t = inter[genes, "t"],
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

rank_age_only <- function(train_counts, train_meta) {
  train_meta$age_group <- factor(train_meta$age_group, levels = c("Young", "Old"))
  design <- model.matrix(~ age_group, data = train_meta)
  if (qr(design)$rank < ncol(design)) stop("Full-data age-only design is not full rank")
  dge <- DGEList(counts = train_counts)
  dge <- normLibSizes(dge, method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tab <- topTable(fit, coef = "age_groupOld", number = Inf, sort.by = "none")
  genes <- rownames(train_counts)
  rownames(tab) <- genes
  r <- data.frame(
    gene = genes,
    age_only_logFC = tab[genes, "logFC"],
    age_only_t = tab[genes, "t"],
    age_only_p_value = tab[genes, "P.Value"],
    age_only_adj_p_value = tab[genes, "adj.P.Val"],
    stringsAsFactors = FALSE
  )
  r$sex_linked_flag <- r$gene %in% sex_linked_genes
  r$module <- ifelse(r$age_only_logFC < 0, "young_high", ifelse(r$age_only_logFC > 0, "old_high", "neutral"))
  r$reliability_pass <- !r$sex_linked_flag & r$module != "neutral" & is.finite(r$age_only_logFC) & is.finite(r$age_only_t)
  r$weight <- abs(r$age_only_logFC) * abs(r$age_only_t)
  r$rank_score <- r$weight * as.numeric(r$reliability_pass)
  r
}

build_sig <- function(ranked, per_direction, selected_genes = NULL, equal_weight = FALSE) {
  r <- ranked[ranked$reliability_pass, , drop = FALSE]
  if (!is.null(selected_genes)) r <- r[r$gene %in% selected_genes, , drop = FALSE]
  y <- r[r$module == "young_high", , drop = FALSE]
  o <- r[r$module == "old_high", , drop = FALSE]
  y <- y[order(-y$rank_score, y$gene), , drop = FALSE]
  o <- o[order(-o$rank_score, o$gene), , drop = FALSE]
  sig <- rbind(head(y, per_direction), head(o, per_direction))
  if (equal_weight && nrow(sig) > 0) sig$weight <- 1
  sig
}

score_signature_expr <- function(signature, expr, metadata_for_calibration) {
  genes <- intersect(signature$gene, rownames(expr))
  sig <- signature[match(genes, signature$gene), , drop = FALSE]
  sig <- sig[is.finite(sig$training_sd) & sig$training_sd > 0, , drop = FALSE]
  genes <- sig$gene
  if (!length(genes)) stop("No usable signature genes for scoring")
  z <- sweep(sweep(expr[genes, , drop = FALSE], 1, sig$training_mean, "-"), 1, sig$training_sd, "/")
  w <- sig$weight
  w[!is.finite(w) | w <= 0] <- 1
  module_score <- function(module) {
    idx <- which(sig$module == module)
    if (!length(idx)) return(rep(NA_real_, ncol(expr)))
    as.numeric(crossprod(w[idx], z[idx, , drop = FALSE]) / sum(w[idx]))
  }
  young_module <- module_score("young_high")
  old_module <- module_score("old_high")
  raw <- young_module - old_module
  young_center <- median(raw[metadata_for_calibration$age_group == "Young"], na.rm = TRUE)
  old_center <- median(raw[metadata_for_calibration$age_group == "Old"], na.rm = TRUE)
  denominator <- young_center - old_center
  calibrated <- (raw - old_center) / denominator
  data.frame(
    mouse = colnames(expr),
    young_module_score = young_module,
    old_module_score = old_module,
    raw_score = raw,
    calibrated_score = calibrated,
    gene_coverage = length(unique(sig$gene)) / nrow(signature),
    weighted_coverage = sum(abs(sig$weight)) / sum(abs(signature$weight)),
    young_reference_center = young_center,
    old_reference_center = old_center,
    calibration_denominator = denominator,
    stringsAsFactors = FALSE
  )
}

freeze_signature <- function(sig, model_name, expr) {
  genes <- sig$gene
  sig$training_mean <- rowMeans(expr[genes, , drop = FALSE])
  sig$training_sd <- apply(expr[genes, , drop = FALSE], 1, sd)
  sig$model <- model_name
  sig$signature_rank <- seq_len(nrow(sig))
  sig[, c("model", "signature_rank", setdiff(names(sig), c("model", "signature_rank"))), drop = FALSE]
}

message("Preparing full-data counts and DE rankings")
prep <- prepare_full(counts, metadata)
filtered_counts <- prep$counts
logcpm <- log2(t(t(filtered_counts) / colSums(filtered_counts)) * 1e6 + 1)
design_full_rank <- qr(model.matrix(~ 0 + age_sex_group, data = metadata))$rank == 4
ranked_factorial <- rank_factorial(filtered_counts, metadata)
ranked_age <- rank_age_only(filtered_counts, metadata)

nested_sigs_path <- file.path(model_comparison_dir, "model_comparison_fold_signatures.csv")
nested_sigs <- read.csv(nested_sigs_path, check.names = FALSE)
primary_nested <- nested_sigs[nested_sigs$model == "factorial_medium_original", , drop = FALSE]
selection_frequency <- sort(table(primary_nested$gene) / length(unique(primary_nested$heldout_mouse)), decreasing = TRUE)
stable_genes_75 <- names(selection_frequency)[selection_frequency >= 0.75]
selection_frequency_df <- data.frame(gene = names(selection_frequency), fold_selection_frequency = as.numeric(selection_frequency), stringsAsFactors = FALSE)
write.csv(selection_frequency_df, file.path(out_dir, "full_data_primary_nested_selection_frequency.csv"), row.names = FALSE)

message("Building frozen signatures")
sig_medium <- freeze_signature(build_sig(ranked_factorial, 50), "factorial_medium_original", logcpm)
sig_large <- freeze_signature(build_sig(ranked_factorial, 100), "factorial_large_original", logcpm)
sig_equal <- freeze_signature(build_sig(ranked_factorial, 50, equal_weight = TRUE), "factorial_medium_equal_weight", logcpm)
sig_stability <- freeze_signature(build_sig(ranked_factorial, 100, selected_genes = stable_genes_75), "factorial_stability_selected", logcpm)
sig_age <- freeze_signature(build_sig(ranked_age, 50), "age_only_de", logcpm)

signatures <- list(
  factorial_medium_original = sig_medium,
  factorial_large_original = sig_large,
  factorial_stability_selected = sig_stability,
  factorial_medium_equal_weight = sig_equal,
  age_only_de = sig_age
)

all_scores <- list()
calibration_rows <- list()
for (model in names(signatures)) {
  sig <- signatures[[model]]
  score <- score_signature_expr(sig, logcpm, metadata)
  score$model <- model
  score <- merge(score, metadata[, c("mouse", "age", "age_months", "age_group", "sex", "n_cells", "pseudobulk_library_size", "pseudobulk_detected_genes")], by = "mouse", all.x = TRUE, sort = FALSE)
  all_scores[[model]] <- score
  calibration_rows[[model]] <- data.frame(
    model = model,
    signature_size = nrow(sig),
    young_high_n = sum(sig$module == "young_high"),
    old_high_n = sum(sig$module == "old_high"),
    young_reference_center = unique(score$young_reference_center),
    old_reference_center = unique(score$old_reference_center),
    calibration_denominator = unique(score$calibration_denominator),
    training_gene_filter_cpm_threshold = cpm_threshold,
    training_gene_filter_min_mice = min_mice_expressed,
    filtered_gene_count = nrow(filtered_counts),
    scoring_transform = "log2(CPM + 1), frozen training gene mean/sd z-score",
    de_normalization = "edgeR TMM + limma voom",
    stringsAsFactors = FALSE
  )
  write.csv(sig, file.path(model_dir, paste0(model, "_signature.csv")), row.names = FALSE)
}
all_signature <- rbind_fill(signatures)
all_scores_df <- do.call(rbind, all_scores)
calibration_df <- do.call(rbind, calibration_rows)

write.csv(all_signature, file.path(model_dir, "facs_v2_full_data_frozen_signatures_all_models.csv"), row.names = FALSE)
write.csv(calibration_df, file.path(model_dir, "facs_v2_full_data_frozen_calibration.csv"), row.names = FALSE)
write.csv(all_scores_df, file.path(out_dir, "full_data_frozen_training_scores.csv"), row.names = FALSE)
saveRDS(list(signatures = signatures, calibration = calibration_df, metadata = metadata, filtered_gene_count = nrow(filtered_counts), design_full_rank = design_full_rank), file.path(model_dir, "facs_v2_full_data_frozen_models.rds"))

model_metrics <- lapply(split(all_scores_df, all_scores_df$model), function(s) {
  data.frame(
    model = unique(s$model),
    auc_young_vs_old = auc_rank(as.integer(s$age_group == "Young"), s$calibrated_score),
    all_age_rho = safe_cor(s$calibrated_score, s$age_months),
    old_only_rho = safe_cor(s$calibrated_score[s$age_group == "Old"], s$age_months[s$age_group == "Old"]),
    young_minus_old_median = median(s$calibrated_score[s$age_group == "Young"], na.rm = TRUE) - median(s$calibrated_score[s$age_group == "Old"], na.rm = TRUE),
    library_rho = safe_cor(s$calibrated_score, s$pseudobulk_library_size),
    detected_rho = safe_cor(s$calibrated_score, s$pseudobulk_detected_genes),
    cell_count_rho = safe_cor(s$calibrated_score, s$n_cells),
    stringsAsFactors = FALSE
  )
})
model_metrics_df <- do.call(rbind, model_metrics)
write.csv(model_metrics_df, file.path(out_dir, "full_data_frozen_training_score_metrics.csv"), row.names = FALSE)

message("Writing parser")
parser_lines <- c(
  "score_facs_v2_youth_model <- function(counts,",
  "                                      signature_csv,",
  "                                      calibration_csv,",
  "                                      model = 'factorial_medium_original',",
  "                                      pseudocount = 1) {",
  "  signature <- read.csv(signature_csv, check.names = FALSE, stringsAsFactors = FALSE)",
  "  calibration <- read.csv(calibration_csv, check.names = FALSE, stringsAsFactors = FALSE)",
  "  signature <- signature[signature$model == model, , drop = FALSE]",
  "  calibration <- calibration[calibration$model == model, , drop = FALSE]",
  "  if (!nrow(signature)) stop('No signature rows for requested model')",
  "  if (!nrow(calibration)) stop('No calibration row for requested model')",
  "  if (is.null(rownames(counts))) stop('counts must have gene rownames')",
  "  if (is.null(colnames(counts))) colnames(counts) <- paste0('sample_', seq_len(ncol(counts)))",
  "  lib <- colSums(counts)",
  "  expr <- log2(t(t(counts) / lib) * 1e6 + pseudocount)",
  "  genes <- intersect(signature$gene, rownames(expr))",
  "  sig <- signature[match(genes, signature$gene), , drop = FALSE]",
  "  sig <- sig[is.finite(sig$training_sd) & sig$training_sd > 0, , drop = FALSE]",
  "  genes <- sig$gene",
  "  if (!length(genes)) stop('No usable signature genes found in counts')",
  "  z <- sweep(sweep(expr[genes, , drop = FALSE], 1, sig$training_mean, '-'), 1, sig$training_sd, '/')",
  "  w <- sig$weight",
  "  w[!is.finite(w) | w <= 0] <- 1",
  "  module_score <- function(module) {",
  "    idx <- which(sig$module == module)",
  "    if (!length(idx)) return(rep(NA_real_, ncol(expr)))",
  "    as.numeric(crossprod(w[idx], z[idx, , drop = FALSE]) / sum(w[idx]))",
  "  }",
  "  young_module <- module_score('young_high')",
  "  old_module <- module_score('old_high')",
  "  raw <- young_module - old_module",
  "  denom <- calibration$calibration_denominator[1]",
  "  calibrated <- (raw - calibration$old_reference_center[1]) / denom",
  "  weighted_coverage <- sum(abs(sig$weight)) / sum(abs(signature$weight))",
  "  data.frame(",
  "    sample_id = colnames(counts),",
  "    model = model,",
  "    young_module_score = young_module,",
  "    old_module_score = old_module,",
  "    raw_score = raw,",
  "    calibrated_score = calibrated,",
  "    gene_coverage = length(unique(sig$gene)) / nrow(signature),",
  "    weighted_coverage = weighted_coverage,",
  "    stringsAsFactors = FALSE",
  "  )",
  "}"
)
parser_path <- file.path(parser_dir, "score_facs_v2_youth_model.R")
writeLines(parser_lines, parser_path)

message("Running parser numerical equivalence test")
source(parser_path)
parser_rows <- list()
for (model in names(signatures)) {
  parsed <- score_facs_v2_youth_model(filtered_counts, file.path(model_dir, "facs_v2_full_data_frozen_signatures_all_models.csv"), file.path(model_dir, "facs_v2_full_data_frozen_calibration.csv"), model = model)
  train <- all_scores_df[all_scores_df$model == model, , drop = FALSE]
  train <- train[match(parsed$sample_id, train$mouse), , drop = FALSE]
  parser_rows[[model]] <- data.frame(
    model = model,
    n_samples = nrow(parsed),
    min_gene_coverage = min(parsed$gene_coverage, na.rm = TRUE),
    min_weighted_coverage = min(parsed$weighted_coverage, na.rm = TRUE),
    max_abs_raw_score_diff = max(abs(parsed$raw_score - train$raw_score), na.rm = TRUE),
    max_abs_calibrated_score_diff = max(abs(parsed$calibrated_score - train$calibrated_score), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
parser_validation_df <- do.call(rbind, parser_rows)
write.csv(parser_validation_df, file.path(out_dir, "parser_numerical_equivalence.csv"), row.names = FALSE)

bootstrap_summary <- if (file.exists(file.path(bootstrap_dir, "donor_bootstrap_ci.csv"))) read.csv(file.path(bootstrap_dir, "donor_bootstrap_ci.csv"), check.names = FALSE) else data.frame()
permutation_summary <- if (file.exists(file.path(permutation_dir, "nested_permutation_medium_summary.csv"))) read.csv(file.path(permutation_dir, "nested_permutation_medium_summary.csv"), check.names = FALSE) else data.frame()

write_model_card <- function() {
  medium_boot <- bootstrap_summary[bootstrap_summary$model == "factorial_medium_original", , drop = FALSE]
  report <- c(
    "# FACS Limb Muscle MSC Youth Score v2 Model Card",
    "",
    "## Model Role",
    "",
    "Primary frozen model: `factorial_medium_original` with 50 young-high and 50 old-high genes. Comparators exported in the same artifact set: `factorial_large_original`, `factorial_stability_selected`, `factorial_medium_equal_weight`, and `age_only_de`.",
    "",
    "This is a FACS TMS Limb Muscle MSC cohort-state score. It should not be described as an externally validated universal MSC aging clock.",
    "",
    "## Training Data",
    "",
    sprintf("- Donors: %d mice", nrow(metadata)),
    sprintf("- Cells represented by pseudobulk: %d", sum(metadata$n_cells)),
    sprintf("- Ages: %s", paste(names(table(metadata$age)), as.integer(table(metadata$age)), sep = "=", collapse = "; ")),
    sprintf("- Sex by age group: %s", paste(names(table(metadata$age_group, metadata$sex)), as.integer(table(metadata$age_group, metadata$sex)), sep = "=", collapse = "; ")),
    sprintf("- Filtered genes: %d", nrow(filtered_counts)),
    sprintf("- Full factorial design full-rank: %s", design_full_rank),
    "",
    "## Frozen Training Procedure",
    "",
    "1. Filter genes at CPM > 1 in at least 2 mice.",
    "2. Use edgeR TMM and limma voom for DE/ranking.",
    "3. Fit full factorial age-by-sex group contrasts and age-only DE.",
    "4. For the primary model, require female/male direction concordance, age-only/common direction concordance, exclusion of listed sex-linked genes, and finite common effects.",
    "5. Rank by `abs(common_logFC) * abs(common_t) * interaction_penalty`.",
    "6. Select top 50 young-high and top 50 old-high genes for Medium Original.",
    "7. Freeze per-gene training mean and sd on full-data log2(CPM + 1).",
    "8. Calibrate raw module difference with full training medians: `(raw - old_reference_center) / (young_reference_center - old_reference_center)`.",
    "",
    "Training-set Young and Old medians are therefore fixed to 1 and 0 by construction; apparent full-data separation is not validation evidence.",
    "",
    "## Validation Evidence Already Completed",
    "",
    "- Donor bootstrap quantified uncertainty for primary and baseline models.",
    "- Formal 999 sex-stratified complete-age-label nested permutation was run for the frozen primary pipeline.",
    "- The primary permutation statistic was `abs_all_age_rho`; AUC distance from 0.5 was supporting only.",
    "",
    "## Key Final Internal Results",
    "",
    paste(capture.output(print(model_metrics_df, row.names = FALSE)), collapse = "\n"),
    "",
    "## Donor Bootstrap: Primary Model",
    "",
    if (nrow(medium_boot)) paste(capture.output(print(medium_boot, row.names = FALSE)), collapse = "\n") else "Bootstrap summary not found.",
    "",
    "## Formal Permutation: Primary Model",
    "",
    if (nrow(permutation_summary)) paste(capture.output(print(permutation_summary, row.names = FALSE)), collapse = "\n") else "Permutation summary not found.",
    "",
    "## Parser Validation",
    "",
    paste(capture.output(print(parser_validation_df, row.names = FALSE)), collapse = "\n"),
    "",
    "## Limitations",
    "",
    "- Strong score-level library-size association remains in the FACS cohort.",
    "- The 3m/18m frozen extension to 24m did not support a monotonic continuous aging trajectory.",
    "- The cohort has only 14 mice and unbalanced 24m sex composition.",
    "- Nested permutation tests label association under the predefined training pipeline, not biological specificity or technical independence.",
    "- No independent external dataset validation has been completed yet.",
    "",
    "## Exported Files",
    "",
    "- `models/facs_v2_full_data_frozen_signatures_all_models.csv`",
    "- `models/facs_v2_full_data_frozen_calibration.csv`",
    "- `models/facs_v2_full_data_frozen_models.rds`",
    "- `parser/score_facs_v2_youth_model.R`",
    "- `full_data_frozen_training_scores.csv`",
    "- `parser_numerical_equivalence.csv`"
  )
  writeLines(report, file.path(out_dir, "FACS_Youth_Score_v2_model_card.md"))
}
write_model_card()

report <- c(
  "# FACS Youth Score v2 Step 18: Full-Data Frozen Model Export",
  "",
  "## Scope",
  "",
  "This step freezes the finalized FACS v2 model artifacts. It does not retune thresholds, change model roles, or inspect cross-assay performance.",
  "",
  "Primary: `factorial_medium_original`.",
  "",
  "Comparators: `factorial_large_original`, `factorial_stability_selected`, `factorial_medium_equal_weight`, `age_only_de`.",
  "",
  "## Model Sizes",
  "",
  paste(capture.output(print(calibration_df[, c("model", "signature_size", "young_high_n", "old_high_n", "calibration_denominator")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Full-Data Apparent Scores",
  "",
  paste(capture.output(print(model_metrics_df, row.names = FALSE)), collapse = "\n"),
  "",
  "These are apparent full-training scores. They are useful for artifact verification and calibration, not validation evidence.",
  "",
  "## Parser Numerical Equivalence",
  "",
  paste(capture.output(print(parser_validation_df, row.names = FALSE)), collapse = "\n"),
  "",
  "## Stability-Selected Definition",
  "",
  sprintf("Stable pool threshold: selected in at least 75%% of Step 14 outer folds. Stable pool size: %d genes. Exported stability-selected signature size: %d genes.", length(stable_genes_75), nrow(sig_stability)),
  "",
  "## Outputs",
  "",
  "- `models/facs_v2_full_data_frozen_signatures_all_models.csv`",
  "- `models/facs_v2_full_data_frozen_calibration.csv`",
  "- `models/facs_v2_full_data_frozen_models.rds`",
  "- `parser/score_facs_v2_youth_model.R`",
  "- `full_data_frozen_training_scores.csv`",
  "- `full_data_frozen_training_score_metrics.csv`",
  "- `parser_numerical_equivalence.csv`",
  "- `FACS_Youth_Score_v2_model_card.md`"
)
writeLines(report, file.path(out_dir, "step18_full_data_frozen_model_export_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step18.txt"))

message("Done")
message("Report: ", file.path(out_dir, "step18_full_data_frozen_model_export_report.md"))
message("Model card: ", file.path(out_dir, "FACS_Youth_Score_v2_model_card.md"))
