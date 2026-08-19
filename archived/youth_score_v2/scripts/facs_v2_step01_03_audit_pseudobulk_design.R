#!/usr/bin/env Rscript

root <- getwd()
local_lib <- file.path(root, ".Rlibs")
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(BPCells)
})

metadata_path <- file.path(root, "data_facs", "limb_muscle_msc", "facs_limb_muscle_msc_young_old_metadata.csv")
matrix_dir <- file.path(root, "data_facs", "limb_muscle_msc", "expression_bpcells_young_old")
out_root <- file.path(root, "outputs", "facs_v2")
qc_dir <- file.path(out_root, "qc")
processed_dir <- file.path(out_root, "processed")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x)]))
  paste(sort(x), collapse = ";")
}

design_diagnostics <- function(design, label) {
  qr_obj <- qr(design)
  singular_values <- svd(design, nu = 0, nv = 0)$d
  condition_number <- if (length(singular_values) == 0 || min(singular_values) == 0) {
    Inf
  } else {
    max(singular_values) / min(singular_values)
  }
  data.frame(
    design = label,
    n_samples = nrow(design),
    n_columns = ncol(design),
    rank = qr_obj$rank,
    full_rank = qr_obj$rank == ncol(design),
    residual_df = nrow(design) - qr_obj$rank,
    condition_number = condition_number,
    columns = paste(colnames(design), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

message("Reading parsed FACS Limb Muscle MSC metadata")
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata$mouse <- metadata$mouse.id
metadata$age_months <- as.integer(sub("m$", "", metadata$age))
metadata$age_group <- factor(metadata$age_group, levels = c("Young", "Old"))
metadata$sex <- factor(metadata$sex, levels = c("female", "male"))
metadata$age_sex_group <- interaction(metadata$sex, metadata$age_group, sep = "_", drop = TRUE)

message("Opening parsed FACS Limb Muscle MSC BPCells matrix")
mat <- open_matrix_dir(matrix_dir)
if (ncol(mat) != nrow(metadata)) {
  stop(sprintf("Matrix/metadata mismatch: %s matrix columns vs %s metadata rows", ncol(mat), nrow(metadata)))
}
if (!identical(as.character(colnames(mat)), as.character(metadata$index))) {
  stop("Matrix column names do not match metadata index order")
}

message("Computing cell-level QC from matrix")
cell_qc <- data.frame(
  index = colnames(mat),
  mouse = metadata$mouse,
  age = metadata$age,
  age_months = metadata$age_months,
  age_group = as.character(metadata$age_group),
  sex = as.character(metadata$sex),
  n_counts_metadata = metadata$n_counts,
  n_genes_metadata = metadata$n_genes,
  n_counts_from_matrix = as.numeric(colSums(mat)),
  n_genes_from_matrix = as.integer(colSums(binarize(mat))),
  stringsAsFactors = FALSE
)
cell_qc$n_counts_metadata_delta <- cell_qc$n_counts_from_matrix - cell_qc$n_counts_metadata
cell_qc$n_genes_metadata_delta <- cell_qc$n_genes_from_matrix - cell_qc$n_genes_metadata
write.csv(cell_qc, file.path(qc_dir, "step01_cell_qc_matrix_metadata_check.csv"), row.names = FALSE)

message("Building mouse-level pseudobulk counts")
mouse_meta <- do.call(rbind, lapply(unique(metadata$mouse), function(mouse_id) {
  idx <- which(metadata$mouse == mouse_id)
  data.frame(
    mouse = mouse_id,
    age = collapse_unique(metadata$age[idx]),
    age_months = unique(metadata$age_months[idx]),
    age_group = collapse_unique(metadata$age_group[idx]),
    sex = collapse_unique(metadata$sex[idx]),
    tissue = collapse_unique(metadata$tissue[idx]),
    cell_ontology_class = collapse_unique(metadata$cell_ontology_class[idx]),
    subtissue = collapse_unique(metadata$subtissue[idx]),
    facs_selection = collapse_unique(metadata$FACS.selection[idx]),
    n_cells = length(idx),
    median_cell_counts = median(cell_qc$n_counts_from_matrix[idx]),
    mean_cell_counts = mean(cell_qc$n_counts_from_matrix[idx]),
    median_cell_detected_genes = median(cell_qc$n_genes_from_matrix[idx]),
    mean_cell_detected_genes = mean(cell_qc$n_genes_from_matrix[idx]),
    stringsAsFactors = FALSE
  )
}))
mouse_meta <- mouse_meta[order(mouse_meta$age_months, mouse_meta$sex, mouse_meta$mouse), ]

pseudobulk <- sapply(mouse_meta$mouse, function(mouse_id) {
  idx <- which(metadata$mouse == mouse_id)
  as.numeric(rowSums(mat[, idx, drop = FALSE]))
})
rownames(pseudobulk) <- rownames(mat)
colnames(pseudobulk) <- mouse_meta$mouse
mouse_meta$pseudobulk_library_size <- as.numeric(colSums(pseudobulk))
mouse_meta$pseudobulk_detected_genes <- as.integer(colSums(pseudobulk > 0))
mouse_meta$age_group <- factor(mouse_meta$age_group, levels = c("Young", "Old"))
mouse_meta$sex <- factor(mouse_meta$sex, levels = c("female", "male"))
mouse_meta$age_sex_group <- interaction(mouse_meta$sex, mouse_meta$age_group, sep = "_", drop = TRUE)

write.csv(mouse_meta, file.path(processed_dir, "facs_v2_limb_msc_mouse_metadata.csv"), row.names = FALSE)
write.csv(
  data.frame(gene = rownames(pseudobulk), pseudobulk, check.names = FALSE),
  file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.csv"),
  row.names = FALSE
)
saveRDS(pseudobulk, file.path(processed_dir, "facs_v2_limb_msc_pseudobulk_counts.rds"))

message("Writing cohort design summaries")
age_summary <- aggregate(
  n_cells ~ age_group + age + age_months,
  data = mouse_meta,
  FUN = sum
)
age_summary$n_mice <- as.integer(aggregate(
  mouse ~ age_group + age + age_months,
  data = mouse_meta,
  FUN = length
)$mouse)
age_summary <- age_summary[order(age_summary$age_months), ]
write.csv(age_summary, file.path(qc_dir, "step01_age_summary.csv"), row.names = FALSE)

age_sex_mouse <- as.data.frame.matrix(table(mouse_meta$age, mouse_meta$sex))
age_sex_mouse$age <- rownames(age_sex_mouse)
age_sex_mouse <- age_sex_mouse[, c("age", setdiff(colnames(age_sex_mouse), "age"))]
write.csv(age_sex_mouse, file.path(qc_dir, "step01_age_by_sex_mouse_counts.csv"), row.names = FALSE)

age_sex_cell <- as.data.frame.matrix(table(metadata$age, metadata$sex))
age_sex_cell$age <- rownames(age_sex_cell)
age_sex_cell <- age_sex_cell[, c("age", setdiff(colnames(age_sex_cell), "age"))]
write.csv(age_sex_cell, file.path(qc_dir, "step01_age_by_sex_cell_counts.csv"), row.names = FALSE)

age_sex_group_counts <- as.data.frame(table(mouse_meta$age_group, mouse_meta$sex))
colnames(age_sex_group_counts) <- c("age_group", "sex", "n_mice")
write.csv(age_sex_group_counts, file.path(qc_dir, "step01_age_group_by_sex_mouse_counts.csv"), row.names = FALSE)

message("Auditing full-data design matrices")
design_age_only <- model.matrix(~ age_group, data = mouse_meta)
design_additive <- model.matrix(~ sex + age_group, data = mouse_meta)
design_factorial <- model.matrix(~ 0 + age_sex_group, data = mouse_meta)
full_design_diagnostics <- rbind(
  design_diagnostics(design_age_only, "age_only_~age_group"),
  design_diagnostics(design_additive, "additive_~sex+age_group"),
  design_diagnostics(design_factorial, "factorial_cell_means_~0+sex_age_group")
)
write.csv(full_design_diagnostics, file.path(qc_dir, "step03_full_design_diagnostics.csv"), row.names = FALSE)
write.csv(design_factorial, file.path(qc_dir, "step03_full_factorial_design_matrix.csv"), row.names = TRUE)

message("Auditing fold-level design support")
fold_rows <- lapply(mouse_meta$mouse, function(heldout) {
  train_meta <- mouse_meta[mouse_meta$mouse != heldout, , drop = FALSE]
  train_meta$age_group <- factor(train_meta$age_group, levels = c("Young", "Old"))
  train_meta$sex <- factor(train_meta$sex, levels = c("female", "male"))
  train_meta$age_sex_group <- interaction(train_meta$sex, train_meta$age_group, sep = "_", drop = TRUE)
  design_fold <- model.matrix(~ 0 + age_sex_group, data = train_meta)
  diag <- design_diagnostics(design_fold, "factorial_cell_means_~0+sex_age_group")
  group_counts <- table(train_meta$age_sex_group)
  expected_groups <- c("female_Young", "male_Young", "female_Old", "male_Old")
  support <- setNames(as.integer(group_counts[expected_groups]), paste0("train_n_", expected_groups))
  support[is.na(support)] <- 0L
  data.frame(
    heldout_mouse = heldout,
    heldout_age = mouse_meta$age[mouse_meta$mouse == heldout],
    heldout_age_months = mouse_meta$age_months[mouse_meta$mouse == heldout],
    heldout_age_group = as.character(mouse_meta$age_group[mouse_meta$mouse == heldout]),
    heldout_sex = as.character(mouse_meta$sex[mouse_meta$mouse == heldout]),
    train_n_mice = nrow(train_meta),
    t(support),
    design_full_rank = diag$full_rank,
    residual_df = diag$residual_df,
    condition_number = diag$condition_number,
    weakly_supported_factorial_cell = any(support < 2),
    stringsAsFactors = FALSE
  )
})
fold_diagnostics <- do.call(rbind, fold_rows)
write.csv(fold_diagnostics, file.path(qc_dir, "step03_fold_design_diagnostics.csv"), row.names = FALSE)

message("Writing audit report")
summary_lines <- c(
  "# FACS v2 Step 01-03: Input Audit, Pseudobulk, and Design Audit",
  "",
  "## Inputs",
  "",
  paste0("- Metadata: `", metadata_path, "`"),
  paste0("- Matrix: `", matrix_dir, "`"),
  "",
  "No files under `data_facs` were modified by this script. All new outputs were written under `outputs/facs_v2`.",
  "",
  "## Input Read Check",
  "",
  paste0("- Matrix dimensions: ", nrow(mat), " genes x ", ncol(mat), " cells"),
  paste0("- Metadata rows: ", nrow(metadata)),
  paste0("- Matrix columns match metadata `index`: ", identical(as.character(colnames(mat)), as.character(metadata$index))),
  paste0("- Metadata `n_counts` exact matches: ", sum(cell_qc$n_counts_metadata_delta == 0), " / ", nrow(cell_qc)),
  paste0("- Metadata `n_genes` exact matches: ", sum(cell_qc$n_genes_metadata_delta == 0), " / ", nrow(cell_qc)),
  "",
  "## Cohort Structure",
  "",
  paste(capture.output(print(mouse_meta[, c("mouse", "age", "age_group", "sex", "n_cells", "pseudobulk_library_size", "pseudobulk_detected_genes")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Age Summary",
  "",
  paste(capture.output(print(age_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Full-Data Design Diagnostics",
  "",
  paste(capture.output(print(full_design_diagnostics, row.names = FALSE)), collapse = "\n"),
  "",
  "## Fold-Level Design Notes",
  "",
  paste0("- Outer LOMO folds audited: ", nrow(fold_diagnostics)),
  paste0("- Full-rank factorial folds: ", sum(fold_diagnostics$design_full_rank), " / ", nrow(fold_diagnostics)),
  paste0("- Folds with any age-by-sex cell supported by <2 training mice: ", sum(fold_diagnostics$weakly_supported_factorial_cell), " / ", nrow(fold_diagnostics)),
  "",
  "Folds holding out one of the two female Young or female Old mice remain algebraically estimable but are weakly supported for sex-specific effects.",
  "",
  "## Outputs",
  "",
  "- `outputs/facs_v2/processed/facs_v2_limb_msc_mouse_metadata.csv`",
  "- `outputs/facs_v2/processed/facs_v2_limb_msc_pseudobulk_counts.csv`",
  "- `outputs/facs_v2/processed/facs_v2_limb_msc_pseudobulk_counts.rds`",
  "- `outputs/facs_v2/qc/step01_cell_qc_matrix_metadata_check.csv`",
  "- `outputs/facs_v2/qc/step01_age_summary.csv`",
  "- `outputs/facs_v2/qc/step01_age_by_sex_mouse_counts.csv`",
  "- `outputs/facs_v2/qc/step01_age_by_sex_cell_counts.csv`",
  "- `outputs/facs_v2/qc/step03_full_design_diagnostics.csv`",
  "- `outputs/facs_v2/qc/step03_fold_design_diagnostics.csv`"
)
writeLines(summary_lines, file.path(qc_dir, "step01_03_audit_pseudobulk_design_report.md"))

message("Done")
message(sprintf("Report: %s", file.path(qc_dir, "step01_03_audit_pseudobulk_design_report.md")))
