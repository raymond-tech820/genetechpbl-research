#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BPCells)
})

matrix_dir <- "data/limb_muscle_msc/expression_bpcells_young_old"
metadata_path <- "data/limb_muscle_msc/limb_muscle_msc_young_old_metadata.csv"
mouse_table_path <- "outputs/qc/mouse_sample_table.csv"
processed_dir <- "data/processed"
qc_dir <- "outputs/qc"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

counts_out <- file.path(processed_dir, "tms_limb_msc_pseudobulk_counts.rds")
metadata_out <- file.path(processed_dir, "tms_limb_msc_pseudobulk_metadata.csv")
counts_csv_out <- file.path(processed_dir, "tms_limb_msc_pseudobulk_counts.csv")
verification_out <- file.path(qc_dir, "step3_pseudobulk_verification.csv")
report_out <- file.path(qc_dir, "step3_pseudobulk_report.md")

message("Opening single-cell BPCells matrix")
single_cell_counts <- open_matrix_dir(matrix_dir)
cell_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
mouse_table <- read.csv(mouse_table_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(as.character(colnames(single_cell_counts)), as.character(cell_metadata$index))) {
  stop("Single-cell matrix colnames do not match metadata index order")
}

cell_metadata$mouse_id <- cell_metadata$mouse.id
mouse_ids <- mouse_table$mouse_id

message("Aggregating raw counts by mouse")
pseudobulk_counts <- matrix(
  0,
  nrow = nrow(single_cell_counts),
  ncol = length(mouse_ids),
  dimnames = list(rownames(single_cell_counts), mouse_ids)
)

for (mouse_id in mouse_ids) {
  cell_idx <- which(cell_metadata$mouse_id == mouse_id)
  if (length(cell_idx) == 0L) {
    stop(sprintf("No cells found for mouse %s", mouse_id))
  }
  pseudobulk_counts[, mouse_id] <- as.numeric(rowSums(single_cell_counts[, cell_idx]))
}

message("Building pseudobulk metadata")
pseudobulk_metadata <- mouse_table[
  ,
  c(
    "mouse_id",
    "age_months",
    "age",
    "age_group",
    "sex",
    "cell_count",
    "total_counts",
    "method",
    "subtissue",
    "library_id"
  )
]
pseudobulk_metadata$sample_id <- pseudobulk_metadata$mouse_id
pseudobulk_metadata <- pseudobulk_metadata[
  ,
  c(
    "sample_id",
    "mouse_id",
    "age_months",
    "age",
    "age_group",
    "sex",
    "cell_count",
    "total_counts",
    "method",
    "subtissue",
    "library_id"
  )
]

message("Running verification checks")
single_cell_gene_totals <- as.numeric(rowSums(single_cell_counts))
pseudobulk_gene_totals <- as.numeric(rowSums(pseudobulk_counts))
gene_total_max_abs_diff <- max(abs(single_cell_gene_totals - pseudobulk_gene_totals))
total_count_single_cell <- sum(single_cell_gene_totals)
total_count_pseudobulk <- sum(pseudobulk_gene_totals)
metadata_total_counts <- mouse_table$total_counts[match(colnames(pseudobulk_counts), mouse_table$mouse_id)]
col_total_max_abs_diff <- max(abs(as.numeric(colSums(pseudobulk_counts)) - metadata_total_counts))

set.seed(20260717)
random_genes <- sample(rownames(pseudobulk_counts), size = 5)
manual_checks <- data.frame(
  gene = character(),
  mouse_id = character(),
  manual_sum = numeric(),
  pseudobulk_count = numeric(),
  matches = logical(),
  stringsAsFactors = FALSE
)
for (gene in random_genes) {
  for (mouse_id in sample(mouse_ids, size = min(3, length(mouse_ids)))) {
    gene_idx <- match(gene, rownames(single_cell_counts))
    cell_idx <- which(cell_metadata$mouse_id == mouse_id)
    manual_sum <- sum(as.numeric(as.matrix(single_cell_counts[gene_idx, cell_idx, drop = FALSE])))
    pb_value <- pseudobulk_counts[gene, mouse_id]
    manual_checks <- rbind(
      manual_checks,
      data.frame(
        gene = gene,
        mouse_id = mouse_id,
        manual_sum = manual_sum,
        pseudobulk_count = pb_value,
        matches = identical(as.numeric(manual_sum), as.numeric(pb_value)),
        stringsAsFactors = FALSE
      )
    )
  }
}
write.csv(manual_checks, file.path(qc_dir, "step3_manual_gene_mouse_checks.csv"), row.names = FALSE)

verification <- data.frame(
  n_genes = nrow(pseudobulk_counts),
  n_mouse_samples = ncol(pseudobulk_counts),
  unique_mice_in_metadata = length(unique(cell_metadata$mouse_id)),
  columns_equal_unique_mice = ncol(pseudobulk_counts) == length(unique(cell_metadata$mouse_id)),
  single_cell_total_counts = total_count_single_cell,
  pseudobulk_total_counts = total_count_pseudobulk,
  total_counts_conserved = total_count_single_cell == total_count_pseudobulk,
  max_abs_gene_total_difference = gene_total_max_abs_diff,
  max_abs_mouse_total_difference_vs_step2 = col_total_max_abs_diff,
  all_manual_gene_mouse_checks_pass = all(manual_checks$matches),
  includes_sex_for_downstream_adjustment = "sex" %in% colnames(pseudobulk_metadata),
  stringsAsFactors = FALSE
)
write.csv(verification, verification_out, row.names = FALSE)

if (!verification$columns_equal_unique_mice) {
  stop("Pseudobulk column count does not match unique mouse count")
}
if (!verification$total_counts_conserved || gene_total_max_abs_diff != 0 || col_total_max_abs_diff != 0) {
  stop("Pseudobulk count conservation check failed")
}
if (!all(manual_checks$matches)) {
  stop("Manual gene-by-mouse checks failed")
}

message("Writing outputs")
saveRDS(pseudobulk_counts, counts_out)
write.csv(pseudobulk_metadata, metadata_out, row.names = FALSE)
counts_csv <- data.frame(gene = rownames(pseudobulk_counts), pseudobulk_counts, check.names = FALSE)
write.csv(counts_csv, counts_csv_out, row.names = FALSE)

report_lines <- c(
  "# Step 3: Mouse-Level Pseudobulk Aggregation",
  "",
  "## Inputs",
  "",
  sprintf("- Single-cell BPCells matrix: `%s`", matrix_dir),
  sprintf("- Cell metadata: `%s`", metadata_path),
  sprintf("- Step 2 mouse table: `%s`", mouse_table_path),
  "",
  "## What This Step Did",
  "",
  "1. Opened the Young/Old Limb_Muscle MSC raw count matrix.",
  "2. Grouped cells by `mouse_id`.",
  "3. Summed raw counts for every gene across all cells from the same mouse.",
  "4. Constructed a gene-by-mouse pseudobulk matrix.",
  "5. Wrote pseudobulk sample metadata with `sex` retained for downstream sex adjustment.",
  "6. Verified count conservation before and after aggregation.",
  "",
  "## Output Dimensions",
  "",
  sprintf("- Genes: %s", nrow(pseudobulk_counts)),
  sprintf("- Mouse-level samples: %s", ncol(pseudobulk_counts)),
  sprintf("- Unique mice in metadata: %s", length(unique(cell_metadata$mouse_id))),
  "",
  "## Verification",
  "",
  sprintf("- Total single-cell counts: %s", total_count_single_cell),
  sprintf("- Total pseudobulk counts: %s", total_count_pseudobulk),
  sprintf("- Max absolute gene-total difference: %s", gene_total_max_abs_diff),
  sprintf("- Max absolute mouse-total difference vs Step 2: %s", col_total_max_abs_diff),
  sprintf("- Manual gene-by-mouse checks pass: %s", all(manual_checks$matches)),
  "",
  "## Sex Adjustment Note",
  "",
  "The pseudobulk metadata keeps `sex` as an explicit column. Downstream model training and differential-expression analyses should include sex adjustment where the design matrix is estimable. Because both young mice are female while old mice include female and male mice, sex remains partially confounded with age and should be reported as a design limitation.",
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", counts_out),
  sprintf("- `%s`", metadata_out),
  sprintf("- `%s`", counts_csv_out),
  sprintf("- `%s`", verification_out),
  "- `outputs/qc/step3_manual_gene_mouse_checks.csv`"
)
writeLines(report_lines, report_out)

message("Step 3 complete")
message(sprintf("Pseudobulk matrix: %s genes x %s mice", nrow(pseudobulk_counts), ncol(pseudobulk_counts)))
message(sprintf("Total counts conserved: %s", total_count_single_cell == total_count_pseudobulk))
