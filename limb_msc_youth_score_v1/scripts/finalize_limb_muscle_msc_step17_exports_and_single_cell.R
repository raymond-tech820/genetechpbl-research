#!/usr/bin/env Rscript

local_lib <- ".Rlibs"
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(BPCells)
  library(ggplot2)
})

signature_versions_path <- "outputs/scores/step12_candidate_signature_versions.csv"
training_scores_path <- "outputs/scores/step12_13_candidate_scores.csv"
training_calibration_path <- "outputs/scores/step13_training_calibration.csv"
step15_summary_path <- "outputs/validation/step15_model_comparison_summary.csv"
step16_summary_path <- "outputs/validation/step16_control_summary.csv"
pseudobulk_logcpm_path <- "data/processed/pseudobulk_logcpm.rds"
cell_matrix_dir <- "data/limb_muscle_msc/expression_bpcells_young_old"
cell_metadata_path <- "data/limb_muscle_msc/limb_muscle_msc_young_old_metadata.csv"
models_dir <- "models"
r_dir <- "R"
single_cell_dir <- "outputs/single_cell"

dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(single_cell_dir, recursive = TRUE, showWarnings = FALSE)

primary_model <- "Medium"
comparator_model <- "Large"
model_version <- "limb_msc_general_youth_score_v1"
minimum_gene_coverage <- 0.8
weighted_coverage_threshold <- 0.8

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x)
  x
}

json_string <- function(x) sprintf('"%s"', json_escape(x))
json_number <- function(x) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) "null" else format(x, digits = 15, scientific = FALSE)
}
json_array_string <- function(x) paste0("[", paste(vapply(x, json_string, character(1)), collapse = ", "), "]")

write_json_object <- function(items, path) {
  lines <- c("{")
  n <- length(items)
  for (i in seq_along(items)) {
    comma <- if (i < n) "," else ""
    lines <- c(lines, sprintf("  %s: %s%s", json_string(names(items)[i]), items[[i]], comma))
  }
  lines <- c(lines, "}")
  writeLines(lines, path)
}

score_signature_matrix <- function(signature, expr_matrix, calibration = NULL, equal_weight = FALSE) {
  sig <- signature[signature$gene %in% rownames(expr_matrix), ]
  if (nrow(sig) == 0) stop("No signature genes found in expression matrix")
  sig <- sig[is.finite(sig$training_sd) & sig$training_sd > 0, ]
  if (equal_weight) sig$weight <- 1
  young_sig <- sig[sig$module == "young_module", ]
  old_sig <- sig[sig$module == "old_module", ]
  score_module <- function(module_sig) {
    z <- sweep(expr_matrix[module_sig$gene, , drop = FALSE], 1, module_sig$training_mean, "-")
    z <- sweep(z, 1, module_sig$training_sd, "/")
    as.numeric(crossprod(abs(module_sig$weight), z) / sum(abs(module_sig$weight)))
  }
  young_score <- score_module(young_sig)
  old_score <- score_module(old_sig)
  raw <- young_score - old_score
  if (is.null(calibration)) {
    calibrated <- rep(NA_real_, length(raw))
    clipped <- rep(NA_real_, length(raw))
  } else {
    calibrated <- (raw - calibration$old_reference_center) / calibration$calibration_denominator
    clipped <- pmin(pmax(calibrated, 0), 1)
  }
  data.frame(
    sample_id = colnames(expr_matrix),
    score_young_module_raw = young_score,
    score_old_module_raw = old_score,
    score_raw = raw,
    youth_score_raw_calibrated = calibrated,
    youth_score_clipped_0_1 = clipped,
    gene_coverage = length(unique(sig$gene)) / nrow(signature),
    weighted_coverage = sum(abs(sig$weight)) / sum(abs(if (equal_weight) rep(1, nrow(signature)) else signature$weight)),
    stringsAsFactors = FALSE
  )
}

message("Reading Step 12-16 outputs")
signature_versions <- read.csv(signature_versions_path, check.names = FALSE, stringsAsFactors = FALSE)
training_scores <- read.csv(training_scores_path, check.names = FALSE, stringsAsFactors = FALSE)
training_calibration <- read.csv(training_calibration_path, check.names = FALSE, stringsAsFactors = FALSE)
step15_summary <- read.csv(step15_summary_path, check.names = FALSE, stringsAsFactors = FALSE)
step16_summary <- read.csv(step16_summary_path, check.names = FALSE, stringsAsFactors = FALSE)

primary_signature <- signature_versions[signature_versions$version == primary_model, ]
large_signature <- signature_versions[signature_versions$version == comparator_model, ]
primary_cal <- training_calibration[training_calibration$version == primary_model, ]
large_cal <- training_calibration[training_calibration$version == comparator_model, ]

final_signature <- data.frame(
  gene = primary_signature$gene,
  direction = primary_signature$direction,
  module = primary_signature$module,
  adjusted_logFC = primary_signature$adjusted_logFC,
  age_rho = primary_signature$age_rho,
  reliability = primary_signature$r_g,
  weight = primary_signature$weight,
  equal_weight = 1,
  training_mean = primary_signature$training_mean,
  training_sd = primary_signature$training_sd,
  FDR = primary_signature$FDR,
  pi_LOMO = primary_signature$pi_LOMO,
  pi_depth = primary_signature$pi_depth,
  pi_sex = primary_signature$pi_sex,
  q_g = primary_signature$q_g,
  primary_model = TRUE,
  stringsAsFactors = FALSE
)
write.csv(final_signature, file.path(models_dir, "limb_msc_general_youth_score_v1_signature.csv"), row.names = FALSE)

equal_weight_signature <- final_signature
equal_weight_signature$weight <- 1
equal_weight_signature$weighting <- "equal_weight_sensitivity"
write.csv(equal_weight_signature, file.path(models_dir, "limb_msc_general_youth_score_v1_signature_equal_weight_medium.csv"), row.names = FALSE)

large_export <- data.frame(
  gene = large_signature$gene,
  direction = large_signature$direction,
  module = large_signature$module,
  adjusted_logFC = large_signature$adjusted_logFC,
  age_rho = large_signature$age_rho,
  reliability = large_signature$r_g,
  weight = large_signature$weight,
  training_mean = large_signature$training_mean,
  training_sd = large_signature$training_sd,
  FDR = large_signature$FDR,
  pi_LOMO = large_signature$pi_LOMO,
  pi_depth = large_signature$pi_depth,
  pi_sex = large_signature$pi_sex,
  q_g = large_signature$q_g,
  primary_model = FALSE,
  comparator_model = TRUE,
  stringsAsFactors = FALSE
)
write.csv(large_export, file.path(models_dir, "limb_msc_general_youth_score_v1_large_comparator_signature.csv"), row.names = FALSE)

metadata_train <- read.csv("data/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv", check.names = FALSE, stringsAsFactors = FALSE)
young_sex_comp <- paste(names(table(metadata_train$sex[metadata_train$age_group == "Young"])),
                        as.integer(table(metadata_train$sex[metadata_train$age_group == "Young"])),
                        sep = ":", collapse = ";")
old_sex_comp <- paste(names(table(metadata_train$sex[metadata_train$age_group == "Old"])),
                      as.integer(table(metadata_train$sex[metadata_train$age_group == "Old"])),
                      sep = ":", collapse = ";")

medium_step15 <- step15_summary[step15_summary$version == primary_model, ]
large_step15 <- step15_summary[step15_summary$version == comparator_model, ]
medium_perm <- step16_summary[step16_summary$version == primary_model &
                                step16_summary$control == "age_label_permutation_training_pipeline", ]
medium_random <- step16_summary[step16_summary$version == primary_model &
                                  step16_summary$control == "expression_matched_random_gene_sets", ]
medium_weight_shuffle <- step16_summary[step16_summary$version == primary_model &
                                          step16_summary$control == "medium_weight_shuffle", ]

calibration_items <- list(
  model_version = json_string(model_version),
  organism = json_string("Mus musculus"),
  tissue = json_string("Limb_Muscle"),
  cell_type = json_string("mesenchymal stem cell"),
  training_dataset = json_string("Tabula Muris Senis droplet Limb_Muscle MSC young_old pseudobulk"),
  normalization_method = json_string("TMM log2 CPM pseudobulk; cell-level export uses module-score sensitivity"),
  primary_model = json_string(primary_model),
  comparator_model = json_string(comparator_model),
  signature_size = json_number(nrow(final_signature)),
  young_high_genes = json_number(sum(final_signature$direction == "young_high")),
  old_high_genes = json_number(sum(final_signature$direction == "old_high")),
  young_reference_center = json_number(primary_cal$young_reference_center),
  old_reference_center = json_number(primary_cal$old_reference_center),
  calibration_denominator = json_number(primary_cal$calibration_denominator),
  minimum_gene_coverage = json_number(minimum_gene_coverage),
  weighted_coverage_threshold = json_number(weighted_coverage_threshold),
  sex_adjustment = json_string("DE and candidate selection used sex-adjusted model ~ sex + age_group; final score is not sex-independent"),
  young_sex_composition = json_string(young_sex_comp),
  old_sex_composition = json_string(old_sex_comp),
  bootstrap_stability = json_string("not_assessed"),
  permutation_scope = json_string("practical_training_pipeline_null"),
  permutation_iterations = json_number(medium_perm$n_successful_permutations),
  random_gene_set_controls = json_number(medium_random$n_controls),
  weight_shuffle_controls = json_number(medium_weight_shuffle$n_controls),
  equal_weight_medium_sensitivity_exported = json_string("yes"),
  primary_lomo_rho_all = json_number(medium_step15$rho_all),
  primary_lomo_rho_old_only = json_number(medium_step15$rho_old_only),
  comparator_lomo_rho_all = json_number(large_step15$rho_all),
  comparator_lomo_rho_old_only = json_number(large_step15$rho_old_only),
  training_mouse_ids = json_array_string(metadata_train$sample_id)
)
write_json_object(calibration_items, file.path(models_dir, "limb_msc_general_youth_score_v1_calibration.json"))

message("Writing R scoring function")
scoring_function <- c(
  "score_limb_msc_youth <- function(expression_logcpm,",
  "                                 signature_path = \"models/limb_msc_general_youth_score_v1_signature.csv\",",
  "                                 calibration_path = \"models/limb_msc_general_youth_score_v1_calibration.json\",",
  "                                 equal_weight = FALSE,",
  "                                 minimum_gene_coverage = NULL,",
  "                                 weighted_coverage_threshold = NULL) {",
  "  parse_json_scalar <- function(path, key) {",
  "    txt <- paste(readLines(path, warn = FALSE), collapse = \"\\n\")",
  "    pattern <- paste0('\"', key, '\"\\\\s*:\\\\s*(\"[^\"]*\"|-?[0-9.]+|null)')",
  "    m <- regexec(pattern, txt)",
  "    hit <- regmatches(txt, m)[[1]]",
  "    if (length(hit) < 2 || hit[2] == \"null\") return(NA)",
  "    val <- hit[2]",
  "    if (startsWith(val, '\"')) return(gsub('^\"|\"$', '', val))",
  "    as.numeric(val)",
  "  }",
  "  sig <- read.csv(signature_path, check.names = FALSE, stringsAsFactors = FALSE)",
  "  if (is.null(rownames(expression_logcpm))) stop('expression_logcpm must have gene rownames')",
  "  if (is.null(colnames(expression_logcpm))) colnames(expression_logcpm) <- paste0('sample_', seq_len(ncol(expression_logcpm)))",
  "  available <- sig$gene %in% rownames(expression_logcpm)",
  "  sig_available <- sig[available, ]",
  "  if (nrow(sig_available) == 0) stop('No signature genes found in expression_logcpm')",
  "  sig_available <- sig_available[is.finite(sig_available$training_sd) & sig_available$training_sd > 0, ]",
  "  if (equal_weight) sig_available$weight <- 1",
  "  if (is.null(minimum_gene_coverage)) minimum_gene_coverage <- as.numeric(parse_json_scalar(calibration_path, 'minimum_gene_coverage'))",
  "  if (is.null(weighted_coverage_threshold)) weighted_coverage_threshold <- as.numeric(parse_json_scalar(calibration_path, 'weighted_coverage_threshold'))",
  "  gene_coverage <- nrow(sig_available) / nrow(sig)",
  "  denominator_weight <- if (equal_weight) nrow(sig) else sum(abs(sig$weight), na.rm = TRUE)",
  "  weighted_coverage <- sum(abs(sig_available$weight), na.rm = TRUE) / denominator_weight",
  "  if (is.finite(minimum_gene_coverage) && gene_coverage < minimum_gene_coverage) warning('Gene coverage below threshold')",
  "  if (is.finite(weighted_coverage_threshold) && weighted_coverage < weighted_coverage_threshold) warning('Weighted coverage below threshold')",
  "  module_score <- function(module_name) {",
  "    module_sig <- sig_available[sig_available$module == module_name, ]",
  "    z <- sweep(expression_logcpm[module_sig$gene, , drop = FALSE], 1, module_sig$training_mean, '-')",
  "    z <- sweep(z, 1, module_sig$training_sd, '/')",
  "    as.numeric(crossprod(abs(module_sig$weight), z) / sum(abs(module_sig$weight)))",
  "  }",
  "  young <- module_score('young_module')",
  "  old <- module_score('old_module')",
  "  raw <- young - old",
  "  old_center <- as.numeric(parse_json_scalar(calibration_path, 'old_reference_center'))",
  "  denom <- as.numeric(parse_json_scalar(calibration_path, 'calibration_denominator'))",
  "  calibrated <- (raw - old_center) / denom",
  "  data.frame(sample_id = colnames(expression_logcpm),",
  "             score_young_module_raw = young,",
  "             score_old_module_raw = old,",
  "             score_raw = raw,",
  "             youth_score_raw_calibrated = calibrated,",
  "             youth_score_clipped_0_1 = pmin(pmax(calibrated, 0), 1),",
  "             gene_coverage = gene_coverage,",
  "             weighted_coverage = weighted_coverage,",
  "             equal_weight = equal_weight,",
  "             stringsAsFactors = FALSE)",
  "}"
)
writeLines(scoring_function, file.path(r_dir, "score_limb_msc_youth.R"))

message("Verifying scoring function against pseudobulk Medium scores")
source(file.path(r_dir, "score_limb_msc_youth.R"))
pseudobulk_logcpm <- readRDS(pseudobulk_logcpm_path)
parser_scores <- score_limb_msc_youth(pseudobulk_logcpm)
expected_medium <- training_scores[training_scores$version == primary_model, ]
verify <- merge(
  expected_medium[, c("sample_id", "score_raw", "youth_score_raw_calibrated")],
  parser_scores[, c("sample_id", "score_raw", "youth_score_raw_calibrated")],
  by = "sample_id",
  suffixes = c("_expected", "_parser")
)
parser_max_abs_raw_diff <- max(abs(verify$score_raw_expected - verify$score_raw_parser))
parser_max_abs_cal_diff <- max(abs(verify$youth_score_raw_calibrated_expected - verify$youth_score_raw_calibrated_parser))

message("Scoring single cells")
cell_metadata <- read.csv(cell_metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
mat <- open_matrix_dir(cell_matrix_dir)
if (!identical(colnames(mat), cell_metadata$index)) {
  stop("Cell metadata index does not match BPCells matrix colnames")
}

sig_genes <- final_signature$gene
missing_cell_genes <- setdiff(sig_genes, rownames(mat))
if (length(missing_cell_genes) > 0) {
  stop(sprintf("Signature genes missing from single-cell matrix: %s", paste(missing_cell_genes, collapse = ",")))
}
sig_counts <- as.matrix(mat[sig_genes, ])
lib_size <- as.numeric(colSums(mat))
names(lib_size) <- colnames(mat)
cell_logcpm <- log2(t(t(sig_counts) / pmax(lib_size[colnames(sig_counts)], 1)) * 1e6 + 1)

score_cell_module <- function(signature, expr, equal_weight = FALSE) {
  sig <- signature[match(rownames(expr), signature$gene), ]
  weights <- if (equal_weight) rep(1, nrow(sig)) else abs(sig$weight)
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0
  young <- sig$module == "young_module"
  old <- sig$module == "old_module"
  young_score <- as.numeric(crossprod(weights[young], z[young, , drop = FALSE]) / sum(weights[young]))
  old_score <- as.numeric(crossprod(weights[old], z[old, , drop = FALSE]) / sum(weights[old]))
  raw <- young_score - old_score
  data.frame(
    index = colnames(expr),
    cell_score_young_module = young_score,
    cell_score_old_module = old_score,
    cell_score_raw = raw,
    equal_weight = equal_weight,
    stringsAsFactors = FALSE
  )
}

cell_scores_weighted <- score_cell_module(final_signature, cell_logcpm, equal_weight = FALSE)
cell_scores_equal <- score_cell_module(final_signature, cell_logcpm, equal_weight = TRUE)
cell_scores <- rbind(
  transform(cell_scores_weighted, score_type = "weighted_medium_module_score"),
  transform(cell_scores_equal, score_type = "equal_weight_medium_module_score")
)
cell_scores <- merge(
  cell_metadata[, c("index", "age", "age_group", "mouse.id", "sex", "n_genes")],
  cell_scores,
  by = "index",
  all.y = TRUE,
  sort = FALSE
)
write.csv(cell_scores, file.path(single_cell_dir, "step17_single_cell_scores.csv"), row.names = FALSE)

aggregate_one <- function(df) {
  aggregate(
    cell_score_raw ~ score_type + mouse.id + age + age_group + sex,
    data = df,
    FUN = function(x) c(n_cells = length(x), mean = mean(x), median = median(x), sd = sd(x))
  )
}
cell_agg <- aggregate_one(cell_scores)
cell_agg <- do.call(data.frame, cell_agg)
names(cell_agg) <- gsub("cell_score_raw\\.", "cell_score_raw_", names(cell_agg))
write.csv(cell_agg, file.path(single_cell_dir, "step17_single_cell_mouse_aggregate_scores.csv"), row.names = FALSE)

pseudobulk_medium <- training_scores[training_scores$version == primary_model, ]
compare_df <- merge(
  cell_agg[cell_agg$score_type == "weighted_medium_module_score", ],
  pseudobulk_medium[, c("mouse_id", "score_raw", "youth_score_raw_calibrated")],
  by.x = "mouse.id",
  by.y = "mouse_id",
  all.x = TRUE
)
comparison_summary <- data.frame(
  score_type = "weighted_medium_module_score",
  n_mice = nrow(compare_df),
  spearman_cell_median_vs_pseudobulk_raw = suppressWarnings(cor(compare_df$cell_score_raw_median, compare_df$score_raw, method = "spearman")),
  pearson_cell_median_vs_pseudobulk_raw = suppressWarnings(cor(compare_df$cell_score_raw_median, compare_df$score_raw, method = "pearson")),
  spearman_cell_mean_vs_pseudobulk_raw = suppressWarnings(cor(compare_df$cell_score_raw_mean, compare_df$score_raw, method = "spearman")),
  bootstrap_stability = "not_assessed",
  cell_score_interpretation = "module_score_sensitivity_not_primary_pseudobulk_parser",
  stringsAsFactors = FALSE
)
write.csv(compare_df, file.path(single_cell_dir, "step17_single_cell_vs_pseudobulk_comparison.csv"), row.names = FALSE)
write.csv(comparison_summary, file.path(single_cell_dir, "step17_single_cell_comparison_summary.csv"), row.names = FALSE)

p <- ggplot(compare_df, aes(x = score_raw, y = cell_score_raw_median, color = age_group, shape = sex)) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, color = "grey35", linewidth = 0.4) +
  scale_color_manual(values = c("Young" = "#D1495B", "Old" = "#4C78A8")) +
  labs(
    title = "Step 17: cell-level aggregate vs pseudobulk Medium score",
    x = "Pseudobulk Medium raw score",
    y = "Per-mouse median single-cell module score",
    color = "Age group",
    shape = "Sex"
  ) +
  theme_classic(base_size = 11)
ggsave(file.path(single_cell_dir, "step17_single_cell_vs_pseudobulk.png"), p, width = 7, height = 5, dpi = 180)

report_lines <- c(
  "# Step 17: Final Medium v1 Export and Single-Cell Sensitivity Score",
  "",
  "## Final Model Decision",
  "",
  "- `primary_model = Medium`",
  "- `comparator_model = Large`",
  "- `bootstrap_stability = not_assessed`",
  "- `permutation_scope = practical_training_pipeline_null`",
  "- Equal-weight Medium sensitivity signature is exported because Step 16 showed that exact weights are not the main performance source.",
  "",
  "## Exports",
  "",
  "- `models/limb_msc_general_youth_score_v1_signature.csv`",
  "- `models/limb_msc_general_youth_score_v1_signature_equal_weight_medium.csv`",
  "- `models/limb_msc_general_youth_score_v1_large_comparator_signature.csv`",
  "- `models/limb_msc_general_youth_score_v1_calibration.json`",
  "- `R/score_limb_msc_youth.R`",
  "",
  "## Parser Verification",
  "",
  sprintf("- Max absolute raw-score difference vs Step 12 Medium: %.6g", parser_max_abs_raw_diff),
  sprintf("- Max absolute calibrated-score difference vs Step 12 Medium: %.6g", parser_max_abs_cal_diff),
  "",
  "## Single-Cell Scoring",
  "",
  "- Single-cell scores use a module-score form on signature-gene log2(CPM+1), z-scored across cells.",
  "- These scores are a sensitivity analysis, not the primary pseudobulk parser score.",
  sprintf("- Cells scored: %s", nrow(cell_metadata)),
  sprintf("- Signature genes used: %s / %s", nrow(final_signature) - length(missing_cell_genes), nrow(final_signature)),
  "",
  "## Single-Cell Aggregate vs Pseudobulk",
  "",
  paste(capture.output(print(comparison_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Outputs",
  "",
  "- `outputs/single_cell/step17_single_cell_scores.csv`",
  "- `outputs/single_cell/step17_single_cell_mouse_aggregate_scores.csv`",
  "- `outputs/single_cell/step17_single_cell_vs_pseudobulk_comparison.csv`",
  "- `outputs/single_cell/step17_single_cell_comparison_summary.csv`",
  "- `outputs/single_cell/step17_single_cell_vs_pseudobulk.png`"
)
writeLines(report_lines, file.path(single_cell_dir, "step17_final_export_and_single_cell_report.md"))

cat("Step 17 complete\n")
cat(sprintf("Parser max raw diff: %.6g\n", parser_max_abs_raw_diff))
cat(sprintf("Parser max calibrated diff: %.6g\n", parser_max_abs_cal_diff))
print(comparison_summary, row.names = FALSE)
