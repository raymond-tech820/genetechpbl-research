#!/usr/bin/env Rscript

set.seed(20260721)

root <- getwd()
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

suppressPackageStartupMessages({
  library(BPCells)
  library(edgeR)
  library(limma)
})

out_root <- file.path(root, "outputs", "facs_v2")
out_dir <- file.path(out_root, "equal_cell_subsampling")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path(root, "data_facs", "limb_muscle_msc", "facs_limb_muscle_msc_young_old_metadata.csv")
matrix_dir <- file.path(root, "data_facs", "limb_muscle_msc", "expression_bpcells_young_old")

n_reps <- 20
medium_per_direction <- 50
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
  labels <- labels[ok]; scores <- scores[ok]
  n_pos <- sum(labels == 1); n_neg <- sum(labels == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

make_mouse_meta <- function(metadata, mice) {
  rows <- lapply(mice, function(mouse_id) {
    idx <- which(metadata$mouse == mouse_id)
    data.frame(
      mouse = mouse_id,
      age = unique(metadata$age[idx]),
      age_months = unique(metadata$age_months[idx]),
      age_group = unique(as.character(metadata$age_group[idx])),
      sex = unique(as.character(metadata$sex[idx])),
      original_n_cells = length(idx),
      stringsAsFactors = FALSE
    )
  })
  meta <- do.call(rbind, rows)
  meta <- meta[order(meta$age_months, meta$sex, meta$mouse), ]
  meta$age_group <- factor(meta$age_group, levels = c("Young", "Old"))
  meta$sex <- factor(meta$sex, levels = c("female", "male"))
  meta$age_sex_group <- interaction(meta$sex, meta$age_group, sep = "_", drop = TRUE)
  meta
}

sample_pseudobulk <- function(mat, metadata, mouse_meta, target_cells) {
  pb <- sapply(mouse_meta$mouse, function(mouse_id) {
    idx <- which(metadata$mouse == mouse_id)
    chosen <- sample(idx, target_cells, replace = FALSE)
    as.numeric(rowSums(mat[, chosen, drop = FALSE]))
  })
  rownames(pb) <- rownames(mat)
  colnames(pb) <- mouse_meta$mouse
  mouse_meta$n_cells <- target_cells
  mouse_meta$pseudobulk_library_size <- as.numeric(colSums(pb))
  mouse_meta$pseudobulk_detected_genes <- as.integer(colSums(pb > 0))
  list(counts = pb, metadata = mouse_meta)
}

rank_genes <- function(train_counts, train_meta) {
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

build_medium <- function(ranked) {
  young <- ranked[ranked$module == "young_high" & ranked$reliability_pass, ]
  old <- ranked[ranked$module == "old_high" & ranked$reliability_pass, ]
  young <- young[order(-young$rank_score, young$gene), ]
  old <- old[order(-old$rank_score, old$gene), ]
  rbind(head(young, medium_per_direction), head(old, medium_per_direction))
}

score_signature <- function(sig, train_expr, held_expr, train_meta) {
  genes <- intersect(sig$gene, rownames(train_expr))
  genes <- intersect(genes, rownames(held_expr))
  sig <- sig[match(genes, sig$gene), ]
  if (length(genes) == 0 || nrow(sig) == 0) return(c(raw = NA, relative = NA))
  mu <- rowMeans(train_expr[genes, , drop = FALSE])
  sdv <- apply(train_expr[genes, , drop = FALSE], 1, sd)
  usable <- is.finite(mu) & is.finite(sdv) & sdv > 0
  genes <- genes[usable]; sig <- sig[usable, ]; mu <- mu[usable]; sdv <- sdv[usable]
  if (length(genes) == 0) return(c(raw = NA, relative = NA))
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
  yc <- median(train_raw[train_meta$age_group == "Young"], na.rm = TRUE)
  oc <- median(train_raw[train_meta$age_group == "Old"], na.rm = TRUE)
  denom <- yc - oc
  rel <- if (is.finite(denom) && abs(denom) > epsilon) (held_raw - oc) / denom else NA_real_
  c(raw = held_raw, relative = rel)
}

run_nested_medium <- function(pb, meta, scenario, replicate_id) {
  score_rows <- list(); sig_rows <- list()
  for (heldout in meta$mouse) {
    held_idx <- which(meta$mouse == heldout)
    train_idx <- setdiff(seq_len(nrow(meta)), held_idx)
    train_meta <- meta[train_idx, ]
    held_meta <- meta[held_idx, ]
    train_meta$age_group <- factor(train_meta$age_group, levels = c("Young", "Old"))
    train_meta$sex <- factor(train_meta$sex, levels = c("female", "male"))
    train_meta$age_sex_group <- factor(interaction(train_meta$sex, train_meta$age_group, sep = "_", drop = TRUE),
      levels = c("female_Young", "male_Young", "female_Old", "male_Old"))
    train_raw <- pb[, train_meta$mouse, drop = FALSE]
    held_raw <- pb[, held_meta$mouse, drop = FALSE]
    keep <- rowSums(t(t(train_raw) / colSums(train_raw)) * 1e6 > cpm_threshold) >= min_mice_expressed
    train_counts <- train_raw[keep, , drop = FALSE]
    held_counts <- held_raw[keep, , drop = FALSE]
    ranked <- rank_genes(train_counts, train_meta)
    sig <- build_medium(ranked)
    sig$scenario <- scenario
    sig$replicate <- replicate_id
    sig$heldout_mouse <- heldout
    sig_rows[[heldout]] <- sig
    deploy_train <- log2(t(t(train_counts) / colSums(train_counts)) * 1e6 + 1)
    deploy_held <- log2(t(t(held_counts) / colSums(held_counts)) * 1e6 + 1)
    sc <- score_signature(sig, deploy_train, deploy_held, train_meta)
    score_rows[[heldout]] <- data.frame(
      scenario = scenario,
      replicate = replicate_id,
      mode = "full_retraining_nested_lomo",
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
      stringsAsFactors = FALSE
    )
  }
  list(scores = do.call(rbind, score_rows), signatures = do.call(rbind, sig_rows))
}

build_frozen_full_data_model <- function(full_counts, full_meta) {
  keep <- rowSums(t(t(full_counts) / colSums(full_counts)) * 1e6 > cpm_threshold) >= min_mice_expressed
  ranked <- rank_genes(full_counts[keep, , drop = FALSE], full_meta)
  sig <- build_medium(ranked)
  expr <- log2(t(t(full_counts[keep, , drop = FALSE]) / colSums(full_counts[keep, , drop = FALSE])) * 1e6 + 1)
  genes <- intersect(sig$gene, rownames(expr))
  sig <- sig[match(genes, sig$gene), ]
  mu <- rowMeans(expr[genes, , drop = FALSE])
  sdv <- apply(expr[genes, , drop = FALSE], 1, sd)
  usable <- is.finite(mu) & is.finite(sdv) & sdv > 0
  sig <- sig[usable, ]; genes <- genes[usable]; mu <- mu[usable]; sdv <- sdv[usable]
  z <- sweep(sweep(expr[genes, , drop = FALSE], 1, mu, "-"), 1, sdv, "/")
  score_one <- function(zv) {
    y <- which(sig$module == "young_high"); o <- which(sig$module == "old_high")
    weighted.mean(zv[y], sig$weight[y], na.rm = TRUE) - weighted.mean(zv[o], sig$weight[o], na.rm = TRUE)
  }
  raw <- vapply(seq_len(ncol(z)), function(j) score_one(z[, j]), numeric(1))
  yc <- median(raw[full_meta$age_group == "Young"], na.rm = TRUE)
  oc <- median(raw[full_meta$age_group == "Old"], na.rm = TRUE)
  list(signature = sig, genes = genes, mu = mu, sd = sdv, young_center = yc, old_center = oc)
}

score_frozen <- function(model, pb, meta, scenario, replicate_id) {
  genes <- intersect(model$genes, rownames(pb))
  expr <- log2(t(t(pb[genes, , drop = FALSE]) / colSums(pb[genes, , drop = FALSE])) * 1e6 + 1)
  sig <- model$signature[match(genes, model$signature$gene), ]
  mu <- model$mu[genes]; sdv <- model$sd[genes]
  z <- sweep(sweep(expr, 1, mu, "-"), 1, sdv, "/")
  y <- which(sig$module == "young_high"); o <- which(sig$module == "old_high")
  raw <- vapply(seq_len(ncol(z)), function(j) {
    weighted.mean(z[y, j], sig$weight[y], na.rm = TRUE) - weighted.mean(z[o, j], sig$weight[o], na.rm = TRUE)
  }, numeric(1))
  rel <- (raw - model$old_center) / (model$young_center - model$old_center)
  data.frame(
    scenario = scenario,
    replicate = replicate_id,
    mode = "frozen_full_data_medium_rescore",
    mouse = meta$mouse,
    age = meta$age,
    age_months = meta$age_months,
    age_group = as.character(meta$age_group),
    sex = as.character(meta$sex),
    score = rel,
    raw_score = raw,
    signature_size = nrow(sig),
    young_high_n = sum(sig$module == "young_high"),
    old_high_n = sum(sig$module == "old_high"),
    library_size = meta$pseudobulk_library_size,
    detected_genes = meta$pseudobulk_detected_genes,
    cell_count = meta$n_cells,
    stringsAsFactors = FALSE
  )
}

summarize_replicates <- function(scores) {
  metric_rows <- lapply(split(scores, paste(scores$scenario, scores$replicate, scores$mode, sep = "::")), function(df) {
    data.frame(
      scenario = df$scenario[1],
      replicate = df$replicate[1],
      mode = df$mode[1],
      auc = auc_rank(as.integer(df$age_group == "Young"), df$score),
      all_age_rho = safe_cor(df$score, df$age_months),
      old_only_rho = safe_cor(df$score[df$age_group == "Old"], df$age_months[df$age_group == "Old"]),
      young_minus_old_median = median(df$score[df$age_group == "Young"], na.rm = TRUE) - median(df$score[df$age_group == "Old"], na.rm = TRUE),
      library_rho = safe_cor(df$score, df$library_size),
      detected_rho = safe_cor(df$score, df$detected_genes),
      cell_count_rho = safe_cor(df$score, df$cell_count),
      zero_signature_folds = sum(df$signature_size == 0),
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, metric_rows)
  summary <- do.call(rbind, lapply(split(metrics, paste(metrics$scenario, metrics$mode, sep = "::")), function(df) {
    numeric_cols <- c("auc", "all_age_rho", "old_only_rho", "young_minus_old_median", "library_rho", "detected_rho", "cell_count_rho", "zero_signature_folds")
    out <- data.frame(scenario = df$scenario[1], mode = df$mode[1], n_replicates = nrow(df))
    for (col in numeric_cols) {
      out[[paste0(col, "_median")]] <- median(df[[col]], na.rm = TRUE)
      out[[paste0(col, "_q025")]] <- quantile(df[[col]], 0.025, na.rm = TRUE)
      out[[paste0(col, "_q975")]] <- quantile(df[[col]], 0.975, na.rm = TRUE)
    }
    out$prop_auc_ge_0_8 <- mean(df$auc >= 0.8, na.rm = TRUE)
    out$prop_all_age_rho_negative <- mean(df$all_age_rho < 0, na.rm = TRUE)
    out$prop_old_only_rho_negative <- mean(df$old_only_rho < 0, na.rm = TRUE)
    out
  }))
  list(metrics = metrics, summary = summary)
}

message("Reading parsed FACS cell matrix and metadata")
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata$mouse <- metadata$mouse.id
metadata$age_months <- as.integer(sub("m$", "", metadata$age))
metadata$age_group <- factor(metadata$age_group, levels = c("Young", "Old"))
metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
mat <- open_matrix_dir(matrix_dir)
stopifnot(identical(as.character(colnames(mat)), as.character(metadata$index)))

all_mice <- make_mouse_meta(metadata, unique(metadata$mouse))
quality_mice <- all_mice$mouse[all_mice$original_n_cells >= 25]
quality_meta_template <- make_mouse_meta(metadata[metadata$mouse %in% quality_mice, ], quality_mice)
scenarios <- list(
  all_mice_14_cells = list(mouse_meta = all_mice, target_cells = 14),
  min25_mice_28_cells = list(mouse_meta = quality_meta_template, target_cells = 28)
)

full_original <- sample_pseudobulk(mat, metadata, all_mice, target_cells = 14)
frozen_model <- build_frozen_full_data_model(full_original$counts, full_original$metadata)

all_score_rows <- list()
all_sig_rows <- list()

for (scenario in names(scenarios)) {
  cfg <- scenarios[[scenario]]
  for (replicate_id in seq_len(n_reps)) {
    message("Scenario ", scenario, " replicate ", replicate_id, "/", n_reps)
    sampled <- sample_pseudobulk(mat, metadata, cfg$mouse_meta, cfg$target_cells)
    nested <- run_nested_medium(sampled$counts, sampled$metadata, scenario, replicate_id)
    frozen <- score_frozen(frozen_model, sampled$counts, sampled$metadata, scenario, replicate_id)
    all_score_rows[[paste(scenario, replicate_id, "nested")]] <- nested$scores
    all_score_rows[[paste(scenario, replicate_id, "frozen")]] <- frozen
    all_sig_rows[[paste(scenario, replicate_id)]] <- nested$signatures
  }
}

scores <- do.call(rbind, all_score_rows)
signatures <- do.call(rbind, all_sig_rows)
summ <- summarize_replicates(scores)

write.csv(scores, file.path(out_dir, "equal_cell_subsampling_scores.csv"), row.names = FALSE)
write.csv(signatures, file.path(out_dir, "equal_cell_subsampling_nested_signatures.csv"), row.names = FALSE)
write.csv(summ$metrics, file.path(out_dir, "equal_cell_subsampling_replicate_metrics.csv"), row.names = FALSE)
write.csv(summ$summary, file.path(out_dir, "equal_cell_subsampling_summary.csv"), row.names = FALSE)

report <- c(
  "# FACS v2 Step 13: Equal-Cell Subsampling Feasibility Audit",
  "",
  "## Scope",
  "",
  "This step follows the negative technical-penalty result. It asks whether donor-level score dependence on library size/cell count persists when each mouse contributes the same number of cells.",
  "",
  paste0("- Replicates per scenario: ", n_reps),
  "- Full retraining mode: every replicate reruns fully nested LOMO Medium signature.",
  "- Frozen rescoring mode: uses a frozen full-data Medium signature and rescales equal-cell pseudobulk.",
  "- No files under `data_facs` were modified.",
  "",
  "## Scenarios",
  "",
  "- `all_mice_14_cells`: all 14 mice, 14 cells per mouse.",
  "- `min25_mice_28_cells`: excludes mice with <25 cells, then uses 28 cells per mouse.",
  "",
  "## Summary",
  "",
  paste(capture.output(print(summ$summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation Guardrails",
  "",
  "This is a 20-replicate feasibility run, not the final 100-500 replicate uncertainty analysis.",
  "If the summary shows persistent high technical correlations or unstable age direction, the cohort-level technical structure remains unresolved.",
  "",
  "## Outputs",
  "",
  "- `equal_cell_subsampling_scores.csv`",
  "- `equal_cell_subsampling_nested_signatures.csv`",
  "- `equal_cell_subsampling_replicate_metrics.csv`",
  "- `equal_cell_subsampling_summary.csv`"
)
writeLines(report, file.path(out_dir, "equal_cell_subsampling_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step13.txt"))

message("Done")
message(sprintf("Report: %s", file.path(out_dir, "equal_cell_subsampling_report.md")))
