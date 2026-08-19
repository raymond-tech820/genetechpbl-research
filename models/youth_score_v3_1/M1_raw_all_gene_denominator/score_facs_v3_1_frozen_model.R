score_facs_v3_1_frozen_model <- function(
  counts,
  signature_csv,
  calibration_csv,
  model = c(
    "M1_raw_all_gene_denominator",
    "M4_raw_all_gene_denominator_pi_0_90"
  ),
  minimum_weighted_coverage = 0.80,
  strict_coverage = TRUE,
  pseudocount = 1
) {
  model <- match.arg(model)
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) {
    stop("counts must have unique gene rownames")
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("counts must contain finite nonnegative values")
  }
  library_size <- colSums(counts)
  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop("every sample must have a positive raw all-gene library size")
  }

  signature <- read.csv(
    signature_csv, check.names = FALSE, stringsAsFactors = FALSE
  )
  calibration <- read.csv(
    calibration_csv, check.names = FALSE, stringsAsFactors = FALSE
  )
  signature <- signature[signature$model == model, , drop = FALSE]
  calibration <- calibration[calibration$model == model, , drop = FALSE]
  if (!nrow(signature)) stop("no signature rows for requested model")
  if (nrow(calibration) != 1L) {
    stop("calibration must contain exactly one row for requested model")
  }
  if (anyDuplicated(signature$gene)) {
    stop("frozen signature contains duplicate genes")
  }

  requested_weight <- sum(signature$weight)
  requested_by_module <- tapply(
    signature$weight, signature$module, sum
  )
  available <- signature$gene %in% rownames(counts) &
    is.finite(signature$training_mean) &
    is.finite(signature$training_sd) & signature$training_sd > 0 &
    is.finite(signature$weight) & signature$weight > 0
  usable_signature <- signature[available, , drop = FALSE]
  if (!nrow(usable_signature)) stop("no usable signature genes in counts")

  available_by_module <- tapply(
    usable_signature$weight, usable_signature$module, sum
  )
  module_names <- c("young_high", "old_high")
  module_weighted_coverage <- setNames(rep(0, 2), module_names)
  for (module in module_names) {
    if (module %in% names(available_by_module) &&
        module %in% names(requested_by_module)) {
      module_weighted_coverage[module] <-
        available_by_module[[module]] / requested_by_module[[module]]
    }
  }
  gene_coverage <- nrow(usable_signature) / nrow(signature)
  weighted_coverage <- sum(usable_signature$weight) / requested_weight
  coverage_pass <- is.finite(weighted_coverage) &&
    weighted_coverage >= minimum_weighted_coverage &&
    all(module_weighted_coverage >= minimum_weighted_coverage)
  if (strict_coverage && !coverage_pass) {
    stop(sprintf(
      paste0(
        "coverage below frozen threshold %.2f: total=%.3f, ",
        "young=%.3f, old=%.3f"
      ),
      minimum_weighted_coverage,
      weighted_coverage,
      module_weighted_coverage[["young_high"]],
      module_weighted_coverage[["old_high"]]
    ))
  }

  expression <- log2(
    t(t(counts) / as.numeric(library_size)) * 1e6 + pseudocount
  )
  genes <- usable_signature$gene
  z <- sweep(
    sweep(
      expression[genes, , drop = FALSE],
      1, usable_signature$training_mean, "-"
    ),
    1, usable_signature$training_sd, "/"
  )
  module_score <- function(module) {
    index <- which(usable_signature$module == module)
    if (!length(index)) return(rep(NA_real_, ncol(counts)))
    as.numeric(
      crossprod(
        usable_signature$weight[index], z[index, , drop = FALSE]
      ) / sum(usable_signature$weight[index])
    )
  }
  young_module <- module_score("young_high")
  old_module <- module_score("old_high")
  raw_score <- young_module - old_module
  denominator <- calibration$calibration_denominator[[1]]
  if (!is.finite(denominator) || denominator <= 0) {
    stop("frozen calibration denominator must be positive")
  }
  calibrated_score <- (
    raw_score - calibration$old_reference_center[[1]]
  ) / denominator

  scores <- data.frame(
    sample = colnames(counts),
    model = model,
    young_module_score = young_module,
    old_module_score = old_module,
    raw_score = raw_score,
    calibrated_score = calibrated_score,
    raw_library_size = as.numeric(library_size),
    gene_coverage = gene_coverage,
    weighted_coverage = weighted_coverage,
    young_high_weighted_coverage =
      module_weighted_coverage[["young_high"]],
    old_high_weighted_coverage =
      module_weighted_coverage[["old_high"]],
    coverage_pass = coverage_pass,
    stringsAsFactors = FALSE
  )
  coverage <- data.frame(
    model = model,
    requested_genes = nrow(signature),
    usable_genes = nrow(usable_signature),
    gene_coverage = gene_coverage,
    weighted_coverage = weighted_coverage,
    young_high_weighted_coverage =
      module_weighted_coverage[["young_high"]],
    old_high_weighted_coverage =
      module_weighted_coverage[["old_high"]],
    minimum_weighted_coverage = minimum_weighted_coverage,
    coverage_pass = coverage_pass,
    stringsAsFactors = FALSE
  )
  list(
    scores = scores,
    coverage = coverage,
    missing_genes = setdiff(signature$gene, usable_signature$gene)
  )
}
