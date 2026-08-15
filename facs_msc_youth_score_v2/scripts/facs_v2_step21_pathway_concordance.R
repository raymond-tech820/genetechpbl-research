#!/usr/bin/env Rscript

# Step 21: pathway-level concordance between FACS v2 Medium and Droplet v1 Medium.
# Uses local Bioconductor GO:BP and Reactome annotation with a unified FACS/Droplet
# filtered-gene universe. No signature/model tuning is performed.

options(stringsAsFactors = FALSE)

root <- "/Users/strangenoah/Desktop/AI+X/Genetech/Youth_score"
out_dir <- file.path(root, "outputs", "facs_v2", "cross_assay_droplet", "pathway_concordance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(GO.db)
  library(reactome.db)
  library(ggplot2)
})

facs_filtered <- readRDS(file.path(root, "outputs", "facs_v2", "processed", "facs_v2_limb_msc_full_data_filtered_counts.rds"))
droplet_filtered <- readRDS(file.path(root, "data_droplet", "processed", "pseudobulk_filtered_counts.rds"))
universe <- sort(intersect(rownames(facs_filtered), rownames(droplet_filtered)))

facs_sig_all <- read.csv(file.path(root, "outputs", "facs_v2", "final_models", "models", "facs_v2_full_data_frozen_signatures_all_models.csv"), check.names = FALSE)
facs_sig <- facs_sig_all[facs_sig_all$model == "factorial_medium_original", , drop = FALSE]
droplet_sig <- read.csv(file.path(root, "models", "limb_msc_general_youth_score_v1_signature.csv"), check.names = FALSE)

gene_sets <- list(
  FACS_v2_young_high = intersect(facs_sig$gene[facs_sig$module == "young_high"], universe),
  FACS_v2_old_high = intersect(facs_sig$gene[facs_sig$module == "old_high"], universe),
  Droplet_v1_young_high = intersect(droplet_sig$gene[droplet_sig$module == "young_module"], universe),
  Droplet_v1_old_high = intersect(droplet_sig$gene[droplet_sig$module == "old_module"], universe)
)

write.csv(data.frame(module = names(gene_sets), n_genes = vapply(gene_sets, length, integer(1))), file.path(out_dir, "pathway_input_module_sizes.csv"), row.names = FALSE)
write.csv(data.frame(gene = universe), file.path(out_dir, "pathway_unified_gene_universe.csv"), row.names = FALSE)

make_go_bp_terms <- function(universe) {
  map <- AnnotationDbi::select(org.Mm.eg.db, keys = universe, keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL"))
  map <- map[!is.na(map$GOALL) & map$ONTOLOGYALL == "BP", c("SYMBOL", "GOALL")]
  map <- unique(map)
  term_info <- AnnotationDbi::select(GO.db, keys = unique(map$GOALL), keytype = "GOID", columns = c("TERM", "ONTOLOGY"))
  term_info <- term_info[term_info$ONTOLOGY == "BP", c("GOID", "TERM")]
  names(term_info) <- c("term_id", "term_name")
  term_genes <- split(map$SYMBOL, map$GOALL)
  data.frame(
    source = "GO:BP",
    term_id = names(term_genes),
    term_name = term_info$term_name[match(names(term_genes), term_info$term_id)],
    genes = vapply(term_genes, function(x) paste(sort(unique(intersect(x, universe))), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
}

make_reactome_terms <- function(universe) {
  sym2entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = universe, keytype = "SYMBOL", columns = "ENTREZID")
  sym2entrez <- unique(sym2entrez[!is.na(sym2entrez$ENTREZID), ])
  entrez_to_sym <- split(sym2entrez$SYMBOL, sym2entrez$ENTREZID)
  ext2path <- as.list(reactomeEXTID2PATHID)
  ext2path <- ext2path[intersect(names(ext2path), names(entrez_to_sym))]
  rows <- list()
  for (eid in names(ext2path)) {
    paths <- ext2path[[eid]]
    syms <- entrez_to_sym[[eid]]
    for (pid in paths) {
      rows[[length(rows) + 1L]] <- data.frame(term_id = as.character(pid), gene = syms, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(source = character(), term_id = character(), term_name = character(), genes = character()))
  m <- unique(do.call(rbind, rows))
  term_genes <- split(m$gene, m$term_id)
  path_names <- as.list(reactomePATHID2NAME)
  data.frame(
    source = "REACTOME",
    term_id = names(term_genes),
    term_name = unname(vapply(names(term_genes), function(id) if (!is.null(path_names[[id]])) path_names[[id]][1] else NA_character_, character(1))),
    genes = vapply(term_genes, function(x) paste(sort(unique(intersect(x, universe))), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
}

ora_one <- function(query_genes, terms, universe, module_name, min_term_size = 5, max_term_size = 500, min_overlap = 1) {
  query_genes <- unique(intersect(query_genes, universe))
  N <- length(universe)
  n <- length(query_genes)
  rows <- vector("list", nrow(terms))
  for (i in seq_len(nrow(terms))) {
    term_genes <- strsplit(terms$genes[i], ";", fixed = TRUE)[[1]]
    term_genes <- unique(term_genes[nzchar(term_genes)])
    term_genes <- intersect(term_genes, universe)
    m <- length(term_genes)
    overlap <- intersect(query_genes, term_genes)
    k <- length(overlap)
    if (m < min_term_size || m > max_term_size || k < min_overlap) next
    p <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
    rows[[i]] <- data.frame(
      module = module_name,
      source = terms$source[i],
      term_id = terms$term_id[i],
      term_name = terms$term_name[i],
      query_size = n,
      universe_size = N,
      term_size = m,
      overlap_n = k,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      p_value = p,
      enrichment_ratio = (k / n) / (m / N),
      ora_score = -log10(pmax(p, .Machine$double.xmin)),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || !nrow(out)) return(data.frame())
  out$adjusted_p_value <- ave(out$p_value, out$source, FUN = function(x) p.adjust(x, method = "BH"))
  out <- out[order(out$adjusted_p_value, out$p_value, -out$overlap_n), ]
  rownames(out) <- NULL
  out
}

message("Building local GO:BP and Reactome term maps")
go_terms <- make_go_bp_terms(universe)
reactome_terms <- make_reactome_terms(universe)
terms <- rbind(go_terms, reactome_terms)
terms$term_genes_n <- lengths(strsplit(terms$genes, ";", fixed = TRUE))
terms <- terms[terms$term_genes_n >= 5 & terms$term_genes_n <= 500 & !is.na(terms$term_name), ]
write.csv(terms[, c("source", "term_id", "term_name", "term_genes_n")], file.path(out_dir, "pathway_term_catalog_used.csv"), row.names = FALSE)

message("Running ORA for four directional modules")
ora_results <- do.call(rbind, lapply(names(gene_sets), function(nm) ora_one(gene_sets[[nm]], terms, universe, nm)))
write.csv(ora_results, file.path(out_dir, "pathway_ora_results_all_modules.csv"), row.names = FALSE)

sig_thr <- 0.10
make_direction_table <- function(ora_results) {
  key <- unique(ora_results[, c("source", "term_id", "term_name")])
  module_names <- names(gene_sets)
  for (m in module_names) {
    x <- ora_results[ora_results$module == m, c("source", "term_id", "adjusted_p_value", "p_value", "enrichment_ratio", "ora_score", "overlap_n", "overlap_genes")]
    names(x)[3:ncol(x)] <- paste(m, names(x)[3:ncol(x)], sep = "__")
    key <- merge(key, x, by = c("source", "term_id"), all.x = TRUE, sort = FALSE)
  }
  key
}

direction_matrix <- make_direction_table(ora_results)
stronger_direction <- function(row, prefix_y, prefix_o) {
  py <- as.numeric(row[[paste0(prefix_y, "__adjusted_p_value")]])
  po <- as.numeric(row[[paste0(prefix_o, "__adjusted_p_value")]])
  if (!is.finite(py)) py <- Inf
  if (!is.finite(po)) po <- Inf
  if (py >= sig_thr && po >= sig_thr) return("not_significant")
  if (py <= po) "young_high" else "old_high"
}
direction_matrix$facs_v2_direction <- apply(direction_matrix, 1, stronger_direction, "FACS_v2_young_high", "FACS_v2_old_high")
direction_matrix$droplet_v1_direction <- apply(direction_matrix, 1, stronger_direction, "Droplet_v1_young_high", "Droplet_v1_old_high")
direction_matrix$concordance <- ifelse(direction_matrix$facs_v2_direction == "not_significant" | direction_matrix$droplet_v1_direction == "not_significant", "not_shared_significant",
  ifelse(direction_matrix$facs_v2_direction == direction_matrix$droplet_v1_direction, paste0("same_", direction_matrix$facs_v2_direction), "opposite_direction"))

signed_score <- function(row, prefix_y, prefix_o) {
  py <- as.numeric(row[[paste0(prefix_y, "__adjusted_p_value")]])
  po <- as.numeric(row[[paste0(prefix_o, "__adjusted_p_value")]])
  sy <- as.numeric(row[[paste0(prefix_y, "__ora_score")]])
  so <- as.numeric(row[[paste0(prefix_o, "__ora_score")]])
  if (!is.finite(py)) py <- Inf
  if (!is.finite(po)) po <- Inf
  if (!is.finite(sy)) sy <- 0
  if (!is.finite(so)) so <- 0
  if (py <= po) sy else -so
}
direction_matrix$facs_v2_signed_ora_score <- apply(direction_matrix, 1, signed_score, "FACS_v2_young_high", "FACS_v2_old_high")
direction_matrix$droplet_v1_signed_ora_score <- apply(direction_matrix, 1, signed_score, "Droplet_v1_young_high", "Droplet_v1_old_high")
write.csv(direction_matrix, file.path(out_dir, "pathway_direction_concordance_matrix.csv"), row.names = FALSE)

shared <- direction_matrix[direction_matrix$concordance %in% c("same_young_high", "same_old_high", "opposite_direction"), ]
shared <- shared[order(shared$concordance, -abs(shared$facs_v2_signed_ora_score) - abs(shared$droplet_v1_signed_ora_score)), ]
write.csv(shared, file.path(out_dir, "pathway_shared_significant_concordance.csv"), row.names = FALSE)

concordance_summary <- data.frame(
  category = names(table(direction_matrix$concordance)),
  n_terms = as.integer(table(direction_matrix$concordance)),
  stringsAsFactors = FALSE
)
write.csv(concordance_summary, file.path(out_dir, "pathway_concordance_summary.csv"), row.names = FALSE)

top_by_module <- do.call(rbind, lapply(split(ora_results, ora_results$module), function(x) head(x[order(x$adjusted_p_value, x$p_value), ], 20)))
write.csv(top_by_module, file.path(out_dir, "pathway_top20_by_module.csv"), row.names = FALSE)

# Plot signed ORA scores for terms significant in at least one model direction and with finite scores.
plot_df <- direction_matrix[direction_matrix$concordance != "not_shared_significant" | abs(direction_matrix$facs_v2_signed_ora_score) > 4 | abs(direction_matrix$droplet_v1_signed_ora_score) > 4, ]
if (nrow(plot_df)) {
  p <- ggplot(plot_df, aes(x = droplet_v1_signed_ora_score, y = facs_v2_signed_ora_score, color = concordance, shape = source)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.2) +
    geom_point(size = 2.2, alpha = 0.8) +
    theme_bw(base_size = 11) +
    labs(title = "Pathway direction concordance", x = "Droplet v1 signed ORA score (+young, -old)", y = "FACS v2 signed ORA score (+young, -old)")
  ggsave(file.path(out_dir, "pathway_signed_ora_concordance.png"), p, width = 7, height = 5, dpi = 200)
}

pkg_versions <- data.frame(
  package = c("AnnotationDbi", "org.Mm.eg.db", "GO.db", "reactome.db"),
  version = c(as.character(packageVersion("AnnotationDbi")), as.character(packageVersion("org.Mm.eg.db")), as.character(packageVersion("GO.db")), as.character(packageVersion("reactome.db"))),
  stringsAsFactors = FALSE
)
write.csv(pkg_versions, file.path(out_dir, "pathway_annotation_package_versions.csv"), row.names = FALSE)

fmt_top <- function(x, n = 8) {
  if (!nrow(x)) return("None")
  paste(capture.output(print(x[seq_len(min(n, nrow(x))), c("module", "source", "term_id", "term_name", "overlap_n", "adjusted_p_value", "enrichment_ratio")], row.names = FALSE)), collapse = "\n")
}

report <- c(
  "# FACS v2 vs Droplet v1 Pathway-Level Concordance",
  "",
  "## Scope",
  "",
  "This analysis asks whether the FACS v2 Medium and Droplet v1 Medium signatures, despite low gene-level overlap, enrich similar directional biological programs on a shared feature background.",
  "",
  "No model parameters were changed. This is interpretive post-hoc annotation of frozen signatures.",
  "",
  "## Gene Set Source",
  "",
  "MSigDB Hallmark via `msigdbr` was attempted but Zenodo downloads repeatedly failed with partial-file errors. g:Profiler API was also attempted but did not connect reliably. To avoid fabricating pathway results, this completed run uses local Bioconductor annotation only:",
  "",
  paste(capture.output(print(pkg_versions, row.names = FALSE)), collapse = "\n"),
  "",
  "Collections used: GO Biological Process and Reactome.",
  "",
  "## Unified Gene Universe",
  "",
  sprintf("FACS filtered genes: %d", nrow(facs_filtered)),
  sprintf("Droplet filtered genes: %d", nrow(droplet_filtered)),
  sprintf("Unified universe intersection: %d", length(universe)),
  "",
  "## Directional Module Sizes",
  "",
  paste(capture.output(print(data.frame(module = names(gene_sets), n_genes = vapply(gene_sets, length, integer(1))), row.names = FALSE)), collapse = "\n"),
  "",
  "## Concordance Summary",
  "",
  paste(capture.output(print(concordance_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Top Shared Significant Terms",
  "",
  if (nrow(shared)) paste(capture.output(print(head(shared[, c("source", "term_id", "term_name", "facs_v2_direction", "droplet_v1_direction", "concordance", "facs_v2_signed_ora_score", "droplet_v1_signed_ora_score")], 20), row.names = FALSE)), collapse = "\n") else "No pathways were significant at FDR < 0.1 in both models under the local GO:BP/Reactome ORA definition.",
  "",
  "## Top Terms By Module",
  "",
  fmt_top(top_by_module, 20),
  "",
  "## Interpretation Guardrails",
  "",
  "- ORA is based on small 50-gene directional modules, so FDR-significant overlap can be sparse.",
  "- Directional concordance is assessed separately for young-high and old-high modules; mixed-direction signatures were not pooled.",
  "- A lack of shared significant pathways does not negate score-level transportability; it means pathway-level mechanism agreement is not strongly supported by this ORA setup.",
  "- Hallmark/SASP/cell-cycle curated collections remain useful future additions if a stable local GMT source is provided.",
  "",
  "## Outputs",
  "",
  "- `pathway_ora_results_all_modules.csv`",
  "- `pathway_direction_concordance_matrix.csv`",
  "- `pathway_shared_significant_concordance.csv`",
  "- `pathway_concordance_summary.csv`",
  "- `pathway_top20_by_module.csv`",
  "- `pathway_term_catalog_used.csv`",
  "- `pathway_signed_ora_concordance.png`"
)
writeLines(report, file.path(out_dir, "pathway_concordance_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_step21_pathway_concordance.txt"))

message("Done")
message("Report: ", file.path(out_dir, "pathway_concordance_report.md"))
