#!/usr/bin/env Rscript

counts_path <- "data/processed/tms_limb_msc_pseudobulk_counts.rds"
metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata.csv"
processed_dir <- "data/processed"
qc_dir <- "outputs/qc"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

filtered_counts_out <- file.path(processed_dir, "pseudobulk_filtered_counts.rds")
filtered_counts_csv_out <- file.path(processed_dir, "pseudobulk_filtered_counts.csv")
gene_qc_out <- file.path(qc_dir, "gene_filter_gene_qc.csv")
summary_out <- file.path(qc_dir, "gene_filter_summary.csv")
report_out <- file.path(qc_dir, "step4_gene_filter_report.md")

cpm_threshold <- 1
min_mice <- 2L
msc_marker_genes <- c("Pdgfra", "Dcn", "Col1a1", "Col1a2", "Col3a1", "Pi16", "Cd34", "Ly6a")

message("Reading pseudobulk counts")
counts <- readRDS(counts_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(colnames(counts), metadata$sample_id)) {
  stop("Pseudobulk count columns do not match pseudobulk metadata sample_id order")
}

library_sizes <- colSums(counts)
if (any(library_sizes <= 0)) {
  stop("One or more pseudobulk samples have non-positive library size")
}

message("Computing CPM and filter metrics")
cpm <- t(t(counts) / library_sizes * 1e6)
detected_mice <- rowSums(cpm > cpm_threshold)
detected_young_mice <- rowSums(cpm[, metadata$age_group == "Young", drop = FALSE] > cpm_threshold)
detected_old_mice <- rowSums(cpm[, metadata$age_group == "Old", drop = FALSE] > cpm_threshold)
total_counts <- rowSums(counts)
mean_cpm <- rowMeans(cpm)
max_cpm <- apply(cpm, 1, max)
keep <- detected_mice >= min_mice

gene_qc <- data.frame(
  gene = rownames(counts),
  total_counts = as.numeric(total_counts),
  mean_cpm = as.numeric(mean_cpm),
  max_cpm = as.numeric(max_cpm),
  mice_with_cpm_gt_1 = as.integer(detected_mice),
  young_mice_with_cpm_gt_1 = as.integer(detected_young_mice),
  old_mice_with_cpm_gt_1 = as.integer(detected_old_mice),
  keep = keep,
  stringsAsFactors = FALSE
)
gene_qc <- gene_qc[order(!gene_qc$keep, -gene_qc$mice_with_cpm_gt_1, -gene_qc$total_counts, gene_qc$gene),]
write.csv(gene_qc, gene_qc_out, row.names = FALSE)

filtered_counts <- counts[keep, , drop = FALSE]
saveRDS(filtered_counts, filtered_counts_out)
filtered_counts_csv <- data.frame(gene = rownames(filtered_counts), filtered_counts, check.names = FALSE)
write.csv(filtered_counts_csv, filtered_counts_csv_out, row.names = FALSE)

marker_status <- gene_qc[match(msc_marker_genes, gene_qc$gene),]
marker_status <- marker_status[!is.na(marker_status$gene),]
write.csv(marker_status, file.path(qc_dir, "gene_filter_msc_marker_status.csv"), row.names = FALSE)

filtered_mean_total_counts <- mean(gene_qc$total_counts[!gene_qc$keep])
retained_mean_total_counts <- mean(gene_qc$total_counts[gene_qc$keep])
coverage_summary <- data.frame(
  cpm_threshold = cpm_threshold,
  min_mice_required = min_mice,
  genes_before_filtering = nrow(counts),
  genes_after_filtering = nrow(filtered_counts),
  genes_removed = nrow(counts) - nrow(filtered_counts),
  fraction_retained = nrow(filtered_counts) / nrow(counts),
  mean_total_counts_retained_genes = retained_mean_total_counts,
  mean_total_counts_removed_genes = filtered_mean_total_counts,
  median_mice_with_cpm_gt_1_retained = median(gene_qc$mice_with_cpm_gt_1[gene_qc$keep]),
  median_mice_with_cpm_gt_1_removed = median(gene_qc$mice_with_cpm_gt_1[!gene_qc$keep]),
  all_samples_retained = ncol(filtered_counts) == nrow(metadata),
  sex_column_available_for_downstream_adjustment = "sex" %in% colnames(metadata),
  stringsAsFactors = FALSE
)
write.csv(coverage_summary, summary_out, row.names = FALSE)

if (nrow(filtered_counts) == 0) {
  stop("Gene filtering removed all genes")
}
if (!identical(colnames(filtered_counts), metadata$sample_id)) {
  stop("Filtered count columns do not match metadata sample_id order")
}

marker_lines <- if (nrow(marker_status) > 0) {
  paste0(
    "- ",
    marker_status$gene,
    ": keep=",
    marker_status$keep,
    ", mice_with_CPM_gt_1=",
    marker_status$mice_with_cpm_gt_1
  )
} else {
  "- None of the marker genes in the audit list were present in the matrix."
}

report_lines <- c(
  "# Step 4: Filter Genes",
  "",
  "## Inputs",
  "",
  sprintf("- Pseudobulk raw counts: `%s`", counts_path),
  sprintf("- Pseudobulk metadata: `%s`", metadata_path),
  "",
  "## What This Step Did",
  "",
  sprintf("1. Computed CPM for every gene in every mouse-level pseudobulk sample."),
  sprintf("2. Counted how many mice had CPM > %s for each gene.", cpm_threshold),
  sprintf("3. Retained genes with CPM > %s in at least %s mice.", cpm_threshold, min_mice),
  "4. Wrote filtered raw counts for downstream normalization and differential-expression analysis.",
  "5. Recorded per-gene coverage metrics and checked a small MSC marker list.",
  "",
  "## Results",
  "",
  sprintf("- Genes before filtering: %s", nrow(counts)),
  sprintf("- Genes after filtering: %s", nrow(filtered_counts)),
  sprintf("- Genes removed: %s", nrow(counts) - nrow(filtered_counts)),
  sprintf("- Fraction retained: %.3f", nrow(filtered_counts) / nrow(counts)),
  sprintf("- Mean total counts among retained genes: %.2f", retained_mean_total_counts),
  sprintf("- Mean total counts among removed genes: %.2f", filtered_mean_total_counts),
  "",
  "## MSC Marker Spot Check",
  "",
  marker_lines,
  "",
  "## Sex Adjustment Note",
  "",
  "Gene filtering is unsupervised with respect to age and sex: the filter uses expression coverage across mice only. The pseudobulk metadata still contains `sex`, and downstream model fitting should include sex adjustment where estimable.",
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", filtered_counts_out),
  sprintf("- `%s`", filtered_counts_csv_out),
  sprintf("- `%s`", gene_qc_out),
  sprintf("- `%s`", summary_out),
  "- `outputs/qc/gene_filter_msc_marker_status.csv`"
)
writeLines(report_lines, report_out)

message("Step 4 complete")
message(sprintf("Genes before filtering: %s", nrow(counts)))
message(sprintf("Genes after filtering: %s", nrow(filtered_counts)))
message(sprintf("Genes removed: %s", nrow(counts) - nrow(filtered_counts)))
