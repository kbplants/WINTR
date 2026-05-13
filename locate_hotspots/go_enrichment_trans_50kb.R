#!/usr/bin/env Rscript
##
## GO enrichment for trans genes per cluster_id from cytoscape_trans_union_50kb.csv.
## Uses clusterProfiler::enricher() with custom GO annotation from 021_locate_hotspot.
## Output: cluster_id, go_enrichment (top GO terms per cluster).
## Run: micromamba run -n r-env Rscript scripts/go_enrichment_trans_50kb.R
##

TOP_N <- 10

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  stop("clusterProfiler not installed. Use: micromamba run -n r-env Rscript ...")
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table", repos = "https://cloud.r-project.org")

library(clusterProfiler)
library(data.table)

WORKDIR <- "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/024_filter_hotspot"
setwd(WORKDIR)

INPUT_CSV  <- file.path("results", "cytoscape_trans_union_50kb.csv")
OUT_CSV    <- file.path("results", "go_enrichment_trans_50kb.csv")
ANNOT_FILE <- "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/021_locate_hotspot/data/Ptrichocarpa_444_v3.1.P14.annotation_info.txt"

to_locus <- function(x) sub("\\.v3\\.1$", "", as.character(x))

message("Reading trans union: ", INPUT_CSV)
dt <- fread(INPUT_CSV)

message("Reading annotation: ", ANNOT_FILE)
ann <- fread(ANNOT_FILE, header = TRUE)
if (substr(names(ann)[1], 1, 1) == "#") names(ann)[1] <- sub("^#", "", names(ann)[1])
locus_col <- "locusName"
go_col    <- "GO"
if (!locus_col %in% names(ann)) stop("Column ", locus_col, " not found.")
if (!go_col %in% names(ann)) stop("Column ", go_col, " not found.")

term2gene_list <- list()
for (i in seq_len(nrow(ann))) {
  locus <- as.character(ann[[locus_col]][i])
  go_str <- as.character(ann[[go_col]][i])
  if (is.na(go_str) || go_str == "" || !grepl("GO:", go_str)) next
  go_ids <- trimws(strsplit(go_str, "[ \t]+")[[1]])
  go_ids <- go_ids[grepl("^GO:[0-9]+", go_ids)]
  for (g in go_ids)
    term2gene_list[[length(term2gene_list) + 1]] <- data.frame(term = g, gene = locus, stringsAsFactors = FALSE)
}
term2gene <- do.call(rbind, term2gene_list)
if (nrow(term2gene) == 0) stop("No GO annotations found in ", ANNOT_FILE)
universe_genes <- unique(term2gene$gene)

term2name <- data.frame(term = unique(term2gene$term), name = unique(term2gene$term), stringsAsFactors = FALSE)
if (requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
  tryCatch({
    go_ids <- unique(term2gene$term)
    names_vec <- vapply(go_ids, function(id) {
      tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[id]]), error = function(e) id)
    }, character(1))
    term2name <- data.frame(term = go_ids, name = names_vec, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

clusters <- unique(dt$cluster_id)
message("Running GO enrichment for ", length(clusters), " clusters (top ", TOP_N, " terms each)...")

out_rows <- list()
for (i in seq_along(clusters)) {
  cid <- clusters[i]
  genes_row <- unique(dt[cluster_id == cid, trans_gene])
  genes_locus <- unique(to_locus(genes_row))
  genes_locus <- intersect(genes_locus, universe_genes)

  if (length(genes_locus) < 5) {
    out_rows[[i]] <- data.frame(cluster_id = cid, go_enrichment = "", stringsAsFactors = FALSE)
    next
  }

  er <- tryCatch(
    enricher(
      gene = genes_locus,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      universe = universe_genes,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 2,
      maxGSSize = 500
    ),
    error = function(e) NULL
  )

  if (is.null(er) || nrow(er@result) == 0) {
    out_rows[[i]] <- data.frame(cluster_id = cid, go_enrichment = "", stringsAsFactors = FALSE)
    next
  }

  res <- as.data.frame(er)
  res <- res[order(res$p.adjust), ]
  n <- min(TOP_N, nrow(res))
  parts <- paste0(res$ID[seq_len(n)], " (", res$Description[seq_len(n)], ")")
  go_str <- paste(parts, collapse = "; ")
  out_rows[[i]] <- data.frame(cluster_id = cid, go_enrichment = go_str, stringsAsFactors = FALSE)

  if (i %% 50 == 0) message("  Cluster ", i, "/", length(clusters))
}

out <- do.call(rbind, out_rows)
write.csv(out, OUT_CSV, row.names = FALSE)
message("Wrote ", OUT_CSV)
message("Done.")
