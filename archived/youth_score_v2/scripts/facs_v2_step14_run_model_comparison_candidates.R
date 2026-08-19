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
out_dir <- file.path(out_root, "model_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.rds"))
metadata <- read.csv(file.path(processed_dir, "facs_v2_limb_msc_mouse_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(colnames(counts), metadata$mouse))
metadata$age_group <- factor(metadata$age_group, levels = c("Young", "Old"))
metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
metadata$age_sex_group <- interaction(metadata$sex, metadata$age_group, sep = "_", drop = TRUE)

cpm_threshold <- 1
min_mice_expressed <- 2
epsilon <- 1e-6
sex_linked_genes <- c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Rps4y1", "Rps4y2", "Sry", "Jarid1d")

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
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

auc_rank <- function(labels, scores) {
  labels <- as.integer(labels)
  ok <- is.finite(scores) & !is.na(labels)
  labels <- labels[ok]; scores <- scores[ok]
  n_pos <- sum(labels == 1); n_neg <- sum(labels == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

prepare_train <- function(count_mat, meta) {
  raw <- count_mat[, meta$mouse, drop = FALSE]
  keep <- rowSums(t(t(raw) / colSums(raw)) * 1e6 > cpm_threshold) >= min_mice_expressed
  train_counts <- raw[keep, , drop = FALSE]
  dge <- DGEList(counts = train_counts)
  dge <- calcNormFactors(dge, method = "TMM")
  list(counts = train_counts, dge = dge, keep = keep)
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

rank_age_only <- function(train_counts, train_meta) {
  dge <- DGEList(counts = train_counts)
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~ age_group, data = train_meta)
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tab <- topTable(fit, coef = "age_groupOld", number = Inf, sort.by = "none")
  genes <- rownames(train_counts)
  rownames(tab) <- genes
  r <- data.frame(
    gene = genes,
    age_only_logFC = tab[genes, "logFC"],
    age_only_t = tab[genes, "t"],
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

score_sig <- function(sig, train_counts, held_counts, train_meta) {
  if (nrow(sig) == 0) return(c(raw = NA, relative = NA, coverage = 0, denom = NA))
  train_expr <- log2(t(t(train_counts) / colSums(train_counts)) * 1e6 + 1)
  held_expr <- log2(t(t(held_counts) / colSums(held_counts)) * 1e6 + 1)
  genes <- intersect(sig$gene, rownames(train_expr))
  sig <- sig[match(genes, sig$gene), , drop = FALSE]
  if (length(genes) == 0 || nrow(sig) == 0) return(c(raw = NA, relative = NA, coverage = 0, denom = NA))
  mu <- rowMeans(train_expr[genes, , drop = FALSE])
  sdv <- apply(train_expr[genes, , drop = FALSE], 1, sd)
  usable <- is.finite(mu) & is.finite(sdv) & sdv > 0
  genes <- genes[usable]; sig <- sig[usable, , drop = FALSE]; mu <- mu[usable]; sdv <- sdv[usable]
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
  c(raw = held_raw, relative = rel, coverage = nrow(sig) / max(1, nrow(sig)), denom = denom)
}

inner_stability_genes <- function(train_meta, count_mat, threshold = 0.75) {
  sigs <- list()
  for (inner_heldout in train_meta$mouse) {
    inner_train <- train_meta[train_meta$mouse != inner_heldout, , drop = FALSE]
    inner_train$age_group <- factor(inner_train$age_group, levels = c("Young", "Old"))
    inner_train$sex <- factor(inner_train$sex, levels = c("female", "male"))
    inner_train$age_sex_group <- factor(interaction(inner_train$sex, inner_train$age_group, sep = "_", drop = TRUE),
      levels = c("female_Young", "male_Young", "female_Old", "male_Old"))
    inner_design <- model.matrix(~ 0 + age_sex_group, data = inner_train)
    if (qr(inner_design)$rank < ncol(inner_design)) {
      next
    }
    prep <- prepare_train(count_mat, inner_train)
    ranked <- rank_factorial(prep$counts, inner_train)
    sig <- build_sig(ranked, 50)
    sigs[[inner_heldout]] <- unique(sig$gene)
  }
  if (length(sigs) == 0) return(character())
  freq <- sort(table(unlist(sigs)), decreasing = TRUE)
  names(freq)[as.numeric(freq) / length(sigs) >= threshold]
}

run_nested_signature_models <- function() {
  score_rows <- list(); sig_rows <- list(); diag_rows <- list()
  models <- c("factorial_medium_original", "factorial_large_original", "factorial_medium_equal_weight", "factorial_stability_selected", "age_only_de")
  for (heldout in metadata$mouse) {
    message("Nested model fold held out: ", heldout)
    held_meta <- metadata[metadata$mouse == heldout, , drop = FALSE]
    train_meta <- metadata[metadata$mouse != heldout, , drop = FALSE]
    train_meta$age_group <- factor(train_meta$age_group, levels = c("Young", "Old"))
    train_meta$sex <- factor(train_meta$sex, levels = c("female", "male"))
    train_meta$age_sex_group <- factor(interaction(train_meta$sex, train_meta$age_group, sep = "_", drop = TRUE),
      levels = c("female_Young", "male_Young", "female_Old", "male_Old"))
    prep <- prepare_train(counts, train_meta)
    held_counts <- counts[prep$keep, heldout, drop = FALSE]
    ranked_factorial <- rank_factorial(prep$counts, train_meta)
    ranked_age <- rank_age_only(prep$counts, train_meta)
    stability_genes <- inner_stability_genes(train_meta, counts, 0.75)
    signatures <- list(
      factorial_medium_original = build_sig(ranked_factorial, 50),
      factorial_large_original = build_sig(ranked_factorial, 100),
      factorial_medium_equal_weight = build_sig(ranked_factorial, 50, equal_weight = TRUE),
      factorial_stability_selected = build_sig(ranked_factorial, 100, selected_genes = stability_genes),
      age_only_de = build_sig(ranked_age, 50)
    )
    for (model in models) {
      sig <- signatures[[model]]
      sig$model <- model; sig$heldout_mouse <- heldout
      sig_rows[[paste(heldout, model)]] <- sig
      sc <- score_sig(sig, prep$counts, held_counts, train_meta)
      score_rows[[paste(heldout, model)]] <- data.frame(
        model = model,
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
        stringsAsFactors = FALSE
      )
    }
    diag_rows[[heldout]] <- data.frame(
      heldout_mouse = heldout,
      retained_genes = sum(prep$keep),
      inner_stability_gene_count = length(stability_genes),
      design_full_rank = qr(model.matrix(~ 0 + age_sex_group, data = train_meta))$rank == 4,
      stringsAsFactors = FALSE
    )
  }
  list(scores = rbind_fill(score_rows), signatures = rbind_fill(sig_rows), diagnostics = rbind_fill(diag_rows))
}

run_pc1_baseline <- function() {
  rows <- list()
  for (heldout in metadata$mouse) {
    train_meta <- metadata[metadata$mouse != heldout, , drop = FALSE]
    held_meta <- metadata[metadata$mouse == heldout, , drop = FALSE]
    prep <- prepare_train(counts, train_meta)
    train_log <- cpm(prep$dge, log = TRUE, prior.count = 1)
    vars <- apply(train_log, 1, var)
    top <- names(sort(vars, decreasing = TRUE))[seq_len(min(2000, length(vars)))]
    pca <- prcomp(t(train_log[top, , drop = FALSE]), center = TRUE, scale. = TRUE)
    held_counts <- counts[prep$keep, heldout, drop = FALSE]
    held_log <- log2(t(t(held_counts) / colSums(held_counts)) * 1e6 + 1)
    held_scaled <- sweep(sweep(t(held_log[top, , drop = FALSE]), 2, pca$center, "-"), 2, pca$scale, "/")
    held_pc1 <- as.numeric(held_scaled %*% pca$rotation[, 1])
    train_pc1 <- pca$x[, 1]
    if (safe_cor(train_pc1, train_meta$age_months) > 0) {
      train_pc1 <- -train_pc1; held_pc1 <- -held_pc1
    }
    yc <- median(train_pc1[train_meta$age_group == "Young"])
    oc <- median(train_pc1[train_meta$age_group == "Old"])
    rows[[heldout]] <- data.frame(
      model = "pc1_baseline",
      mouse = heldout,
      age = held_meta$age,
      age_months = held_meta$age_months,
      age_group = as.character(held_meta$age_group),
      sex = as.character(held_meta$sex),
      score = (held_pc1 - oc) / (yc - oc),
      raw_score = held_pc1,
      signature_size = NA_integer_,
      young_high_n = NA_integer_,
      old_high_n = NA_integer_,
      library_size = held_meta$pseudobulk_library_size,
      detected_genes = held_meta$pseudobulk_detected_genes,
      cell_count = held_meta$n_cells,
      calibration_denominator = yc - oc,
      gene_coverage = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_models <- function(scores, signatures = NULL) {
  rows <- list()
  for (model in unique(scores$model)) {
    s <- scores[scores$model == model, , drop = FALSE]
    sig <- if (!is.null(signatures) && "model" %in% colnames(signatures)) signatures[signatures$model == model, , drop = FALSE] else data.frame()
    sets <- if (nrow(sig)) split(sig$gene, sig$heldout_mouse) else list()
    jac <- NA_real_; genes75 <- NA_integer_
    if (length(sets) > 1) {
      pairs <- combn(names(sets), 2, simplify = FALSE)
      jac <- median(vapply(pairs, function(p) length(intersect(sets[[p[1]]], sets[[p[2]]])) / length(union(sets[[p[1]]], sets[[p[2]]])), numeric(1)), na.rm = TRUE)
      freq <- table(unlist(sets))
      genes75 <- sum(as.numeric(freq) / length(sets) >= 0.75)
    }
    rows[[model]] <- data.frame(
      model = model,
      auc = auc_rank(as.integer(s$age_group == "Young"), s$score),
      all_age_rho = safe_cor(s$score, s$age_months),
      old_only_rho = safe_cor(s$score[s$age_group == "Old"], s$age_months[s$age_group == "Old"]),
      young_minus_old_median = median(s$score[s$age_group == "Young"], na.rm = TRUE) - median(s$score[s$age_group == "Old"], na.rm = TRUE),
      library_rho = safe_cor(s$score, s$library_size),
      detected_rho = safe_cor(s$score, s$detected_genes),
      cell_count_rho = safe_cor(s$score, s$cell_count),
      median_signature_size = median(s$signature_size, na.rm = TRUE),
      median_jaccard = jac,
      genes_selected_ge_75pct = genes75,
      zero_signature_folds = sum(s$signature_size == 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

message("Running nested model-comparison candidates")
nested <- run_nested_signature_models()
pc1 <- run_pc1_baseline()
scores <- rbind(nested$scores, pc1)

elastic_path <- file.path(out_root, "baselines", "elastic_net", "facs_v2_elastic_net_summary.csv")
elastic_available <- file.exists(elastic_path)

summary <- summarize_models(scores, nested$signatures)
summary$role <- c(
  factorial_medium_original = "Primary predefined",
  factorial_large_original = "Large comparator",
  factorial_medium_equal_weight = "Weight sensitivity",
  factorial_stability_selected = "Stability comparator",
  age_only_de = "Biological baseline",
  pc1_baseline = "Technical-axis baseline"
)[summary$model]

if (elastic_available) {
  elastic <- read.csv(elastic_path, stringsAsFactors = FALSE)
} else {
  elastic <- data.frame(note = "Elastic-net baseline summary not found in expected facs_v2 path")
}

write.csv(scores, file.path(out_dir, "model_comparison_nested_scores.csv"), row.names = FALSE)
write.csv(nested$signatures, file.path(out_dir, "model_comparison_fold_signatures.csv"), row.names = FALSE)
write.csv(nested$diagnostics, file.path(out_dir, "model_comparison_fold_diagnostics.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "model_comparison_summary.csv"), row.names = FALSE)
write.csv(elastic, file.path(out_dir, "elastic_net_reference_summary.csv"), row.names = FALSE)

report <- c(
  "# FACS v2 Step 14: Main and Comparator Model Runs",
  "",
  "## Scope",
  "",
  "Technical-tightening is closed. This step runs predefined main/comparator models without retuning thresholds, signature sizes, penalty strengths, or donor exclusions.",
  "",
  "Models run here:",
  "",
  "- Factorial Medium Original",
  "- Factorial Large Original",
  "- Factorial Medium Equal Weight",
  "- Factorial Stability-Selected using inner-fold >=75% selection inside each outer fold",
  "- Age-Only Signed-DE baseline",
  "- Nested PC1 baseline",
  "",
  "Kaile frozen model was not run because no Kaile model file is currently available. GSE176206 was not used.",
  "",
  "No files under `data_facs` were modified.",
  "",
  "## Summary",
  "",
  paste(capture.output(print(summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Notes",
  "",
  "- The PC1 model is a technical-axis comparator, not a biological Youth Score.",
  "- The stability-selected comparator uses inner stability within each outer fold to avoid held-out leakage.",
  "- Model roles remain predefined; this table is for credibility assessment, not post-hoc primary-model replacement.",
  "",
  "## Outputs",
  "",
  "- `model_comparison_nested_scores.csv`",
  "- `model_comparison_fold_signatures.csv`",
  "- `model_comparison_fold_diagnostics.csv`",
  "- `model_comparison_summary.csv`",
  "- `elastic_net_reference_summary.csv`"
)
writeLines(report, file.path(out_dir, "model_comparison_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step14.txt"))

message("Done")
message(sprintf("Report: %s", file.path(out_dir, "model_comparison_report.md")))
