score_facs_v2_1_youth_model <- function(counts,
                                      signature_csv,
                                      calibration_csv,
                                      model = 'factorial_medium_original',
                                      pseudocount = 1) {
  signature <- read.csv(signature_csv, check.names = FALSE, stringsAsFactors = FALSE)
  calibration <- read.csv(calibration_csv, check.names = FALSE, stringsAsFactors = FALSE)
  signature <- signature[signature$model == model, , drop = FALSE]
  calibration <- calibration[calibration$model == model, , drop = FALSE]
  if (!nrow(signature)) stop('No signature rows for requested model')
  if (!nrow(calibration)) stop('No calibration row for requested model')
  if (is.null(rownames(counts))) stop('counts must have gene rownames')
  if (is.null(colnames(counts))) colnames(counts) <- paste0('sample_', seq_len(ncol(counts)))
  lib <- colSums(counts)
  expr <- log2(t(t(counts) / lib) * 1e6 + pseudocount)
  genes <- intersect(signature$gene, rownames(expr))
  sig <- signature[match(genes, signature$gene), , drop = FALSE]
  sig <- sig[is.finite(sig$training_sd) & sig$training_sd > 0, , drop = FALSE]
  genes <- sig$gene
  if (!length(genes)) stop('No usable signature genes found in counts')
  z <- sweep(sweep(expr[genes, , drop = FALSE], 1, sig$training_mean, '-'), 1, sig$training_sd, '/')
  w <- sig$weight
  w[!is.finite(w) | w <= 0] <- 1
  module_score <- function(module) {
    idx <- which(sig$module == module)
    if (!length(idx)) return(rep(NA_real_, ncol(expr)))
    as.numeric(crossprod(w[idx], z[idx, , drop = FALSE]) / sum(w[idx]))
  }
  young_module <- module_score('young_high')
  old_module <- module_score('old_high')
  raw <- young_module - old_module
  denom <- calibration$calibration_denominator[1]
  calibrated <- (raw - calibration$old_reference_center[1]) / denom
  weighted_coverage <- sum(abs(sig$weight)) / sum(abs(signature$weight))
  data.frame(
    sample_id = colnames(counts),
    model = model,
    young_module_score = young_module,
    old_module_score = old_module,
    raw_score = raw,
    calibrated_score = calibrated,
    gene_coverage = length(unique(sig$gene)) / nrow(signature),
    weighted_coverage = weighted_coverage,
    stringsAsFactors = FALSE
  )
}
