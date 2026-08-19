#!/usr/bin/env Rscript

# FACS Youth Score v2 Step 17: formal practical nested permutation null.
# Primary scope: Factorial Medium Original only.
# Formal definition:
# - Sex-stratified permutation of the complete 3m/18m/24m age label.
# - Each permutation fully reruns outer LOMO feature selection, TMM, voom,
#   factorial DE, ranking, signature selection, weights, calibration, and scoring.
# - Primary statistic: abs(all-age Spearman rho).
# - Supporting statistic: abs(AUC - 0.5).
# - Purpose: test age-label association beyond randomized labels, not technical independence.

options(stringsAsFactors = FALSE)

root <- "/Users/strangenoah/Desktop/AI+X/Genetech/Youth_score"
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

n_perm <- as.integer(Sys.getenv("FACS_V2_N_PERM", "999"))
base_seed <- as.integer(Sys.getenv("FACS_V2_PERM_BASE_SEED", "20260721"))
out_root <- file.path(root, "outputs", "facs_v2")
processed_dir <- file.path(out_root, "processed")
out_dir <- file.path(out_root, "permutation", "formal_medium_999")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.rds"))
metadata0 <- read.csv(file.path(processed_dir, "facs_v2_limb_msc_mouse_metadata.csv"), check.names = FALSE)
stopifnot(identical(colnames(counts), metadata0$mouse))

metadata0$age_group <- ifelse(metadata0$age_months == 3, "Young", "Old")
metadata0$age_group <- factor(metadata0$age_group, levels = c("Young", "Old"))
metadata0$sex <- factor(metadata0$sex, levels = c("female", "male"))

cpm_threshold <- 1
min_mice_expressed <- 2
epsilon <- 1e-6
sex_linked_genes <- c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Rps4y1", "Rps4y2", "Sry", "Jarid1d")

score_file <- file.path(out_dir, "nested_permutation_medium_scores.csv")
diag_file <- file.path(out_dir, "nested_permutation_medium_fold_diagnostics.csv")
metric_file <- file.path(out_dir, "nested_permutation_medium_metrics.csv")
summary_file <- file.path(out_dir, "nested_permutation_medium_summary.csv")
report_file <- file.path(out_dir, "nested_permutation_medium_report.md")

append_csv <- function(df, path) {
  write.table(df, path, sep = ",", row.names = FALSE, col.names = !file.exists(path), append = file.exists(path), quote = TRUE)
}

md5_text <- function(x) {
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(x, tf, useBytes = TRUE)
  unname(tools::md5sum(tf))
}

assignment_id <- function(meta) {
  z <- meta[order(meta$mouse), c("mouse", "sex", "age", "age_months", "age_group")]
  paste(paste(z$mouse, z$sex, z$age, z$age_months, as.character(z$age_group), sep = ":"), collapse = "|")
}

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

prepare_train <- function(count_mat, meta) {
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
  if (qr(design_age)$rank < ncol(design_age) || qr(design_fact)$rank < ncol(design_fact)) return(NULL)

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
    female_logFC = female[genes, "logFC"],
    male_logFC = male[genes, "logFC"],
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

build_sig <- function(ranked, per_direction = 50) {
  if (is.null(ranked)) return(ranked)
  r <- ranked[ranked$reliability_pass, , drop = FALSE]
  y <- r[r$module == "young_high", , drop = FALSE]
  o <- r[r$module == "old_high", , drop = FALSE]
  y <- y[order(-y$rank_score, y$gene), , drop = FALSE]
  o <- o[order(-o$rank_score, o$gene), , drop = FALSE]
  rbind(head(y, per_direction), head(o, per_direction))
}

score_sig <- function(sig, train_counts, held_counts, train_meta) {
  if (is.null(sig) || nrow(sig) == 0) return(c(raw = NA, relative = NA, denom = NA))
  train_expr <- log2(t(t(train_counts) / colSums(train_counts)) * 1e6 + 1)
  held_expr <- log2(t(t(held_counts) / colSums(held_counts)) * 1e6 + 1)
  genes <- intersect(sig$gene, rownames(train_expr))
  sig <- sig[match(genes, sig$gene), , drop = FALSE]
  if (!length(genes) || nrow(sig) == 0) return(c(raw = NA, relative = NA, denom = NA))
  mu <- rowMeans(train_expr[genes, , drop = FALSE])
  sdv <- apply(train_expr[genes, , drop = FALSE], 1, sd)
  usable <- is.finite(mu) & is.finite(sdv) & sdv > 0
  genes <- genes[usable]
  sig <- sig[usable, , drop = FALSE]
  mu <- mu[usable]
  sdv <- sdv[usable]
  if (!length(genes) || nrow(sig) == 0) return(c(raw = NA, relative = NA, denom = NA))
  ztrain <- sweep(sweep(train_expr[genes, , drop = FALSE], 1, mu, "-"), 1, sdv, "/")
  zheld <- (as.numeric(held_expr[genes, 1]) - mu) / sdv
  w <- sig$weight
  w[!is.finite(w) | w <= 0] <- 1
  mod_score <- function(z, module) {
    idx <- which(sig$module == module)
    if (!length(idx)) return(NA_real_)
    weighted.mean(z[idx], w[idx], na.rm = TRUE)
  }
  train_raw <- vapply(seq_len(ncol(ztrain)), function(j) mod_score(ztrain[, j], "young_high") - mod_score(ztrain[, j], "old_high"), numeric(1))
  held_raw <- mod_score(zheld, "young_high") - mod_score(zheld, "old_high")
  yc <- median(train_raw[train_meta$age_group == "Young"], na.rm = TRUE)
  oc <- median(train_raw[train_meta$age_group == "Old"], na.rm = TRUE)
  denom <- yc - oc
  rel <- if (is.finite(denom) && abs(denom) > epsilon) (held_raw - oc) / denom else NA_real_
  c(raw = held_raw, relative = rel, denom = denom)
}

permute_age_within_sex <- function(meta, seed) {
  set.seed(seed)
  out <- meta
  out$observed_age <- out$age
  out$observed_age_months <- out$age_months
  out$observed_age_group <- as.character(out$age_group)
  for (sx in levels(out$sex)) {
    idx <- which(out$sex == sx)
    perm_idx <- sample(idx, length(idx), replace = FALSE)
    out$age[idx] <- meta$age[perm_idx]
    out$age_months[idx] <- meta$age_months[perm_idx]
  }
  out$age_group <- factor(ifelse(out$age_months == 3, "Young", "Old"), levels = c("Young", "Old"))
  out
}

run_nested_medium <- function(meta, permutation_id, seed, label_assignment_hash, label_assignment_id) {
  meta$age_group <- factor(meta$age_group, levels = c("Young", "Old"))
  meta$sex <- factor(meta$sex, levels = c("female", "male"))
  rows <- vector("list", nrow(meta))
  diags <- vector("list", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    heldout <- meta$mouse[i]
    held_meta <- meta[i, , drop = FALSE]
    train_meta <- meta[meta$mouse != heldout, , drop = FALSE]
    train_meta$age_group <- factor(train_meta$age_group, levels = c("Young", "Old"))
    train_meta$sex <- factor(train_meta$sex, levels = c("female", "male"))
    train_meta$age_sex_group <- factor(interaction(train_meta$sex, train_meta$age_group, sep = "_", drop = TRUE), levels = c("female_Young", "male_Young", "female_Old", "male_Old"))
    design <- model.matrix(~ 0 + age_sex_group, data = train_meta)
    design_full_rank <- qr(design)$rank == ncol(design)
    if (!design_full_rank) {
      rows[[i]] <- data.frame(permutation_id, seed, label_assignment_hash, mouse = heldout, age = held_meta$age, age_months = held_meta$age_months,
        age_group = as.character(held_meta$age_group), sex = as.character(held_meta$sex), score = NA_real_, raw_score = NA_real_, signature_size = NA_integer_,
        young_high_n = NA_integer_, old_high_n = NA_integer_, library_size = held_meta$pseudobulk_library_size,
        detected_genes = held_meta$pseudobulk_detected_genes, cell_count = held_meta$n_cells, calibration_denominator = NA_real_)
      diags[[i]] <- data.frame(permutation_id, seed, label_assignment_hash, heldout_mouse = heldout, design_full_rank = FALSE, retained_genes = NA_integer_, signature_size = NA_integer_, weak_support_fold = NA)
      next
    }
    prep <- prepare_train(counts, train_meta)
    ranked <- rank_factorial(prep$counts, train_meta)
    sig <- build_sig(ranked, 50)
    sc <- score_sig(sig, prep$counts, counts[prep$keep, heldout, drop = FALSE], train_meta)
    ct <- table(train_meta$age_sex_group)
    rows[[i]] <- data.frame(
      permutation_id, seed, label_assignment_hash,
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
      calibration_denominator = sc[["denom"]]
    )
    diags[[i]] <- data.frame(permutation_id, seed, label_assignment_hash, heldout_mouse = heldout, design_full_rank = TRUE, retained_genes = sum(prep$keep), signature_size = nrow(sig), weak_support_fold = any(ct <= 1))
  }
  scores <- do.call(rbind, rows)
  diagnostics <- do.call(rbind, diags)
  attr(scores, "label_assignment_id") <- label_assignment_id
  list(scores = scores, diagnostics = diagnostics)
}

calc_metrics <- function(scores, diagnostics, permutation_id, seed, type, label_assignment_id, label_assignment_hash, valid, failure_reason) {
  ok_scores <- scores[is.finite(scores$score), , drop = FALSE]
  auc <- auc_rank(as.integer(ok_scores$age_group == "Young"), ok_scores$score)
  rho <- safe_cor(ok_scores$score, ok_scores$age_months)
  old <- ok_scores[ok_scores$age_group == "Old", , drop = FALSE]
  data.frame(
    permutation_id = permutation_id,
    seed = seed,
    label_assignment_id = label_assignment_id,
    label_assignment_hash = label_assignment_hash,
    type = type,
    valid = valid,
    failure_reason = failure_reason,
    n_scored_folds = nrow(ok_scores),
    auc = auc,
    auc_distance_from_half = ifelse(is.finite(auc), abs(auc - 0.5), NA_real_),
    all_age_rho = rho,
    abs_all_age_rho = ifelse(is.finite(rho), abs(rho), NA_real_),
    old_only_rho = safe_cor(old$score, old$age_months),
    young_minus_old_median = median(ok_scores$score[ok_scores$age_group == "Young"], na.rm = TRUE) - median(ok_scores$score[ok_scores$age_group == "Old"], na.rm = TRUE),
    library_rho = safe_cor(ok_scores$score, ok_scores$library_size),
    detected_rho = safe_cor(ok_scores$score, ok_scores$detected_genes),
    cell_count_rho = safe_cor(ok_scores$score, ok_scores$cell_count),
    zero_signature_folds = ifelse(nrow(diagnostics), sum(diagnostics$signature_size == 0, na.rm = TRUE), NA_integer_),
    weak_support_folds = ifelse(nrow(diagnostics), sum(diagnostics$weak_support_fold, na.rm = TRUE), NA_integer_),
    rank_deficient_folds = ifelse(nrow(diagnostics), sum(!diagnostics$design_full_rank, na.rm = TRUE), NA_integer_),
    na_score_folds = ifelse(nrow(scores), sum(!is.finite(scores$score)), NA_integer_)
  )
}

empty_scores <- function(permutation_id, seed, hash) {
  data.frame(permutation_id = integer(), seed = integer(), label_assignment_hash = character(), mouse = character(), age = character(), age_months = numeric(), age_group = character(), sex = character(), score = numeric(), raw_score = numeric(), signature_size = integer(), young_high_n = integer(), old_high_n = integer(), library_size = numeric(), detected_genes = numeric(), cell_count = integer(), calibration_denominator = numeric())
}
empty_diags <- function(permutation_id, seed, hash) {
  data.frame(permutation_id = integer(), seed = integer(), label_assignment_hash = character(), heldout_mouse = character(), design_full_rank = logical(), retained_genes = integer(), signature_size = integer(), weak_support_fold = logical())
}

completed <- if (file.exists(metric_file)) read.csv(metric_file, check.names = FALSE) else data.frame()
completed_ids <- if (nrow(completed)) completed$permutation_id else integer()
used_hashes <- if (nrow(completed)) completed$label_assignment_hash else character()
observed_assignment_id <- assignment_id(metadata0)
observed_hash <- md5_text(observed_assignment_id)

if (!(0L %in% completed_ids)) {
  message("Running observed nested primary model")
  res <- tryCatch(run_nested_medium(metadata0, 0L, base_seed, observed_hash, observed_assignment_id), error = function(e) e)
  if (inherits(res, "error")) {
    scores <- empty_scores(0L, base_seed, observed_hash)
    diags <- empty_diags(0L, base_seed, observed_hash)
    met <- calc_metrics(scores, diags, 0L, base_seed, "observed", observed_assignment_id, observed_hash, FALSE, conditionMessage(res))
  } else {
    scores <- res$scores
    diags <- res$diagnostics
    met <- calc_metrics(scores, diags, 0L, base_seed, "observed", observed_assignment_id, observed_hash, TRUE, "")
  }
  append_csv(scores, score_file)
  append_csv(diags, diag_file)
  append_csv(met, metric_file)
  completed_ids <- c(completed_ids, 0L)
  used_hashes <- c(used_hashes, observed_hash)
}

for (p in seq_len(n_perm)) {
  if (p %in% completed_ids) next
  seed <- base_seed + p * 1009L
  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    perm_meta <- permute_age_within_sex(metadata0, seed + attempt - 1L)
    aid <- assignment_id(perm_meta)
    ahash <- md5_text(aid)
    if (!(ahash %in% used_hashes)) break
    if (attempt > 100L) stop("Could not generate a unique sex-stratified assignment for permutation ", p)
  }
  if (p %% 25L == 0L || p == 1L) message("Permutation ", p, "/", n_perm, " seed=", seed, " assignment=", substr(ahash, 1, 10))
  res <- tryCatch(run_nested_medium(perm_meta, p, seed, ahash, aid), error = function(e) e)
  if (inherits(res, "error")) {
    scores <- empty_scores(p, seed, ahash)
    diags <- empty_diags(p, seed, ahash)
    met <- calc_metrics(scores, diags, p, seed, "sex_stratified_complete_age_label_permutation", aid, ahash, FALSE, conditionMessage(res))
  } else {
    scores <- res$scores
    diags <- res$diagnostics
    met <- calc_metrics(scores, diags, p, seed, "sex_stratified_complete_age_label_permutation", aid, ahash, TRUE, "")
  }
  append_csv(scores, score_file)
  append_csv(diags, diag_file)
  append_csv(met, metric_file)
  used_hashes <- c(used_hashes, ahash)
}

metrics_df <- read.csv(metric_file, check.names = FALSE)
obs <- metrics_df[metrics_df$permutation_id == 0 & metrics_df$valid, , drop = FALSE]
null <- metrics_df[metrics_df$permutation_id != 0 & metrics_df$valid, , drop = FALSE]

emp_upper <- function(col, obs_val) {
  vals <- null[[col]][is.finite(null[[col]])]
  if (!length(vals) || !is.finite(obs_val)) return(c(p = NA_real_, n = length(vals)))
  c(p = (1 + sum(vals >= obs_val)) / (length(vals) + 1), n = length(vals))
}

summary_rows <- list(
  data.frame(metric = "abs_all_age_rho", role = "primary", observed = obs$abs_all_age_rho, null_median = median(null$abs_all_age_rho, na.rm = TRUE), null_ci_low = quantile(null$abs_all_age_rho, 0.025, names = FALSE, na.rm = TRUE), null_ci_high = quantile(null$abs_all_age_rho, 0.975, names = FALSE, na.rm = TRUE), empirical_p = emp_upper("abs_all_age_rho", obs$abs_all_age_rho)["p"], p_mode = "upper_tail_abs_spearman", valid_null_n = emp_upper("abs_all_age_rho", obs$abs_all_age_rho)["n"]),
  data.frame(metric = "auc_distance_from_half", role = "supporting", observed = obs$auc_distance_from_half, null_median = median(null$auc_distance_from_half, na.rm = TRUE), null_ci_low = quantile(null$auc_distance_from_half, 0.025, names = FALSE, na.rm = TRUE), null_ci_high = quantile(null$auc_distance_from_half, 0.975, names = FALSE, na.rm = TRUE), empirical_p = emp_upper("auc_distance_from_half", obs$auc_distance_from_half)["p"], p_mode = "upper_tail_abs_auc_minus_0.5", valid_null_n = emp_upper("auc_distance_from_half", obs$auc_distance_from_half)["n"]),
  data.frame(metric = "young_minus_old_median", role = "recorded_not_primary", observed = obs$young_minus_old_median, null_median = median(null$young_minus_old_median, na.rm = TRUE), null_ci_low = quantile(null$young_minus_old_median, 0.025, names = FALSE, na.rm = TRUE), null_ci_high = quantile(null$young_minus_old_median, 0.975, names = FALSE, na.rm = TRUE), empirical_p = emp_upper("young_minus_old_median", obs$young_minus_old_median)["p"], p_mode = "upper_tail_positive_difference", valid_null_n = emp_upper("young_minus_old_median", obs$young_minus_old_median)["n"])
)
perm_summary <- do.call(rbind, summary_rows)
write.csv(perm_summary, summary_file, row.names = FALSE)

fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = 3, format = "f"))
report <- c(
  "# FACS Youth Score v2 Step 17: Formal Nested Permutation Null for Primary Medium Model",
  "",
  sprintf("Attempted permutations requested: %d", n_perm),
  sprintf("Valid null permutations completed: %d", nrow(null)),
  sprintf("Failed permutations: %d", sum(metrics_df$permutation_id != 0 & !metrics_df$valid)),
  "",
  "## Locked Definition",
  "",
  "- Model: Factorial Medium Original only.",
  "- Permutation: complete 3m/18m/24m age labels are reassigned within sex strata, then `age_months` and `age_group` are regenerated from the permuted full age label.",
  "- Female label composition is preserved as 2 x 3m and 2 x 18m.",
  "- Male label composition is preserved as 4 x 3m, 2 x 18m, and 4 x 24m.",
  "- Every permutation reruns the full nested LOMO training-side workflow, including gene filtering, TMM, voom, factorial DE, ranking, signature selection, weights, calibration, and held-out scoring.",
  "- Primary statistic: `abs_all_age_rho = abs(Spearman(score, age_months))`.",
  "- Supporting statistic: `auc_distance_from_half = abs(AUC - 0.5)`.",
  "- Empirical p-value: `(1 + number of null statistics >= observed statistic) / (valid_null_n + 1)`.",
  "",
  "This is a practical training-pipeline permutation null. It tests whether the predefined nested factorial pipeline produces stronger age-label association than sex-stratified randomized age labels. It does not test technical independence, external validity, or whether the score is a continuous biological aging clock.",
  "",
  "## Summary",
  "",
  "| Metric | Role | Observed | Null median | Null 95% interval | Empirical p | Valid null n |",
  "|---|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(perm_summary))) {
  z <- perm_summary[i, ]
  report <- c(report, sprintf("| %s | %s | %s | %s | [%s, %s] | %s | %d |", z$metric, z$role, fmt(z$observed), fmt(z$null_median), fmt(z$null_ci_low), fmt(z$null_ci_high), fmt(z$empirical_p), as.integer(z$valid_null_n)))
}
report <- c(report, "", "## Diagnostics", "", sprintf("- Total metric rows: %d", nrow(metrics_df)), sprintf("- Duplicate assignment hashes among valid null permutations: %d", sum(duplicated(null$label_assignment_hash))), sprintf("- Observed zero-signature folds: %s", obs$zero_signature_folds), sprintf("- Observed rank-deficient folds: %s", obs$rank_deficient_folds), "", "## Outputs", "", "- `nested_permutation_medium_scores.csv`", "- `nested_permutation_medium_fold_diagnostics.csv`", "- `nested_permutation_medium_metrics.csv`", "- `nested_permutation_medium_summary.csv`")
writeLines(report, report_file)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step17_permutation.txt"))

message("Done")
message("Report: ", report_file)
