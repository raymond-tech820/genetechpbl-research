#!/usr/bin/env Rscript

metadata_path <- "data/processed/tms_limb_msc_pseudobulk_metadata.csv"
norm_factors_path <- "outputs/qc/tmm_normalization_factors.csv"
processed_dir <- "data/processed"
qc_dir <- "outputs/qc"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

labels_out <- file.path(processed_dir, "tms_limb_msc_training_labels.csv")
metadata_labeled_out <- file.path(processed_dir, "tms_limb_msc_pseudobulk_metadata_labeled.csv")
verification_out <- file.path(qc_dir, "step7_label_verification.csv")
report_out <- file.path(qc_dir, "step7_age_label_report.md")

message("Reading pseudobulk metadata")
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
norm_factors <- read.csv(norm_factors_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!identical(metadata$sample_id, norm_factors$sample_id)) {
  stop("Pseudobulk metadata order does not match normalization-factor table")
}

valid_young <- metadata$age_months == 3 & metadata$age_group == "Young"
valid_old <- metadata$age_months %in% c(18, 21, 24) & metadata$age_group == "Old"
if (!all(valid_young | valid_old)) {
  stop("Unexpected age/age_group combination found")
}

metadata$young_label <- ifelse(metadata$age_months == 3, 1L, 0L)
metadata$old_label <- 1L - metadata$young_label
metadata$age_group_model <- factor(metadata$age_group, levels = c("Old", "Young"))
metadata$sex_model <- factor(metadata$sex)
metadata$age_months_continuous <- metadata$age_months
metadata$age_months_scaled <- as.numeric(scale(metadata$age_months))
metadata$effective_library_size <- norm_factors$effective_library_size
metadata$raw_library_size <- norm_factors$raw_library_size
metadata$tmm_norm_factor <- norm_factors$tmm_norm_factor

design <- model.matrix(~ sex_model + age_group_model, data = metadata)
design_rank <- qr(design)$rank
design_full_rank <- design_rank == ncol(design)

labels <- metadata[
  ,
  c(
    "sample_id",
    "mouse_id",
    "age_months",
    "age_months_continuous",
    "age_months_scaled",
    "age",
    "age_group",
    "young_label",
    "old_label",
    "sex",
    "cell_count",
    "raw_library_size",
    "effective_library_size",
    "tmm_norm_factor",
    "library_id"
  )
]
write.csv(labels, labels_out, row.names = FALSE)
write.csv(metadata, metadata_labeled_out, row.names = FALSE)

verification <- data.frame(
  n_samples = nrow(metadata),
  young_label_1_count = sum(metadata$young_label == 1),
  young_label_0_count = sum(metadata$young_label == 0),
  young_mice = sum(metadata$age_group == "Young"),
  old_mice = sum(metadata$age_group == "Old"),
  young_ages = paste(sort(unique(metadata$age_months[metadata$young_label == 1])), collapse = ";"),
  old_ages = paste(sort(unique(metadata$age_months[metadata$young_label == 0])), collapse = ";"),
  all_3m_are_young_label_1 = all(metadata$young_label[metadata$age_months == 3] == 1),
  all_18_21_24m_are_young_label_0 = all(metadata$young_label[metadata$age_months %in% c(18, 21, 24)] == 0),
  sex_adjusted_design_formula = "~ sex + age_group",
  sex_adjusted_design_columns = paste(colnames(design), collapse = ";"),
  sex_adjusted_design_rank = design_rank,
  sex_adjusted_design_ncol = ncol(design),
  sex_adjusted_design_full_rank = design_full_rank,
  age_sex_interaction_used = FALSE,
  stringsAsFactors = FALSE
)
write.csv(verification, verification_out, row.names = FALSE)

if (!design_full_rank) {
  stop("Sex-adjusted age-group design matrix is not full rank")
}

age_table <- capture.output(print(table(metadata$age_group, metadata$age_months)))
label_table <- capture.output(print(table(metadata$age_group, metadata$young_label)))
sex_table <- capture.output(print(table(metadata$age_group, metadata$sex)))

report_lines <- c(
  "# Step 7: Define Age Labels",
  "",
  "## Inputs",
  "",
  sprintf("- Pseudobulk metadata: `%s`", metadata_path),
  sprintf("- TMM normalization factors: `%s`", norm_factors_path),
  "",
  "## What This Step Did",
  "",
  "1. Defined the binary young label for the main comparison.",
  "2. Assigned `young_label = 1` to 3-month mice.",
  "3. Assigned `young_label = 0` to 18-, 21-, and 24-month mice.",
  "4. Retained chronological age in months as `age_months_continuous` for monotonic trend checks.",
  "5. Added scaled chronological age as `age_months_scaled` for later diagnostics if needed.",
  "6. Retained `sex`, cell count, raw library size, and TMM effective library size for downstream adjustment and diagnostics.",
  "7. Checked that the future sex-adjusted design `~ sex + age_group` is full rank.",
  "",
  "## Label Definition",
  "",
  "```text",
  "young_label = 1 for 3m",
  "young_label = 0 for 18m, 21m, 24m",
  "age_group reference for DE = Old",
  "age_group contrast for DE = Young vs Old",
  "```",
  "",
  "## Verification",
  "",
  sprintf("- Total mouse-level samples: %s", nrow(metadata)),
  sprintf("- Young-label samples: %s", sum(metadata$young_label == 1)),
  sprintf("- Old-label samples: %s", sum(metadata$young_label == 0)),
  sprintf("- Young ages represented: %s", verification$young_ages),
  sprintf("- Old ages represented: %s", verification$old_ages),
  sprintf("- Sex-adjusted design full rank: %s", design_full_rank),
  sprintf("- Design columns: %s", paste(colnames(design), collapse = ", ")),
  "",
  "Age by age-month table:",
  "",
  "```text",
  age_table,
  "```",
  "",
  "Age group by young_label table:",
  "",
  "```text",
  label_table,
  "```",
  "",
  "Age group by sex table:",
  "",
  "```text",
  sex_table,
  "```",
  "",
  "## Sex Adjustment Note",
  "",
  "This step does not fit a model, but it prepares labels for the guidance-compliant downstream model `~ sex + age_group`. No age-by-sex interaction is included or planned because young male mice are unavailable.",
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", labels_out),
  sprintf("- `%s`", metadata_labeled_out),
  sprintf("- `%s`", verification_out)
)
writeLines(report_lines, report_out)

message("Step 7 complete")
message(sprintf("Samples: %s; young_label=1: %s; young_label=0: %s", nrow(metadata), sum(metadata$young_label == 1), sum(metadata$young_label == 0)))
message(sprintf("Sex-adjusted design full rank: %s", design_full_rank))
