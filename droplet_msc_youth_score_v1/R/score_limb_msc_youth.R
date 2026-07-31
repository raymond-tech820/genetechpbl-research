score_limb_msc_youth <- function(expression_logcpm,
                                 signature_path = "models/limb_msc_general_youth_score_v1_signature.csv",
                                 calibration_path = "models/limb_msc_general_youth_score_v1_calibration.json",
                                 equal_weight = FALSE,
                                 minimum_gene_coverage = NULL,
                                 weighted_coverage_threshold = NULL) {
  parse_json_scalar <- function(path, key) {
    txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
    pattern <- paste0('"', key, '"\\s*:\\s*("[^"]*"|-?[0-9.]+|null)')
    m <- regexec(pattern, txt)
    hit <- regmatches(txt, m)[[1]]
    if (length(hit) < 2 || hit[2] == "null") return(NA)
    val <- hit[2]
    if (startsWith(val, '"')) return(gsub('^"|"$', '', val))
    as.numeric(val)
  }
  sig <- read.csv(signature_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(rownames(expression_logcpm))) stop('expression_logcpm must have gene rownames')
  if (is.null(colnames(expression_logcpm))) colnames(expression_logcpm) <- paste0('sample_', seq_len(ncol(expression_logcpm)))
  available <- sig$gene %in% rownames(expression_logcpm)
  sig_available <- sig[available, ]
  if (nrow(sig_available) == 0) stop('No signature genes found in expression_logcpm')
  sig_available <- sig_available[is.finite(sig_available$training_sd) & sig_available$training_sd > 0, ]
  if (equal_weight) sig_available$weight <- 1
  if (is.null(minimum_gene_coverage)) minimum_gene_coverage <- as.numeric(parse_json_scalar(calibration_path, 'minimum_gene_coverage'))
  if (is.null(weighted_coverage_threshold)) weighted_coverage_threshold <- as.numeric(parse_json_scalar(calibration_path, 'weighted_coverage_threshold'))
  gene_coverage <- nrow(sig_available) / nrow(sig)
  denominator_weight <- if (equal_weight) nrow(sig) else sum(abs(sig$weight), na.rm = TRUE)
  weighted_coverage <- sum(abs(sig_available$weight), na.rm = TRUE) / denominator_weight
  if (is.finite(minimum_gene_coverage) && gene_coverage < minimum_gene_coverage) warning('Gene coverage below threshold')
  if (is.finite(weighted_coverage_threshold) && weighted_coverage < weighted_coverage_threshold) warning('Weighted coverage below threshold')
  module_score <- function(module_name) {
    module_sig <- sig_available[sig_available$module == module_name, ]
    z <- sweep(expression_logcpm[module_sig$gene, , drop = FALSE], 1, module_sig$training_mean, '-')
    z <- sweep(z, 1, module_sig$training_sd, '/')
    as.numeric(crossprod(abs(module_sig$weight), z) / sum(abs(module_sig$weight)))
  }
  young <- module_score('young_module')
  old <- module_score('old_module')
  raw <- young - old
  old_center <- as.numeric(parse_json_scalar(calibration_path, 'old_reference_center'))
  denom <- as.numeric(parse_json_scalar(calibration_path, 'calibration_denominator'))
  calibrated <- (raw - old_center) / denom
  data.frame(sample_id = colnames(expression_logcpm),
             score_young_module_raw = young,
             score_old_module_raw = old,
             score_raw = raw,
             youth_score_raw_calibrated = calibrated,
             youth_score_clipped_0_1 = pmin(pmax(calibrated, 0), 1),
             gene_coverage = gene_coverage,
             weighted_coverage = weighted_coverage,
             equal_weight = equal_weight,
             stringsAsFactors = FALSE)
}
