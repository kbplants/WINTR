#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# Usage:
#   Rscript scripts/build_triplets_20kb_winter.R .
# or:
#   Rscript scripts/build_triplets_20kb_winter.R /path/to/032_cis_mediation_fdr20
#
# Inputs (in data/):
#   - trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt
#   - Ptrichocarpa_444_v3.1.gene_exons.gff3
# Output (in results/winter/):
#   - cis_trans_triplets_all_20kb_winter.tsv
#
# For each SNP in the winter trans-eQTL (≥50 genes, FDR ≤ 0.2) file, this script
# finds genes with TSS within ±10 kb (using gene coordinates in the GFF3 file),
# then builds SNP–window-gene–trans triplets.

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "."

trans_file <- file.path(base_dir, "data", "trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt")
gff_file <- file.path(base_dir, "data", "Ptrichocarpa_444_v3.1.gene_exons.gff3")
out_file <- file.path(base_dir, "results", "winter", "cis_trans_triplets_all_20kb_winter.tsv")

message("Loading winter trans-eQTL results (≥50 genes, FDR ≤ 0.2)...")
trans_dt <- fread(trans_file)

if (!all(c("snps", "gene") %in% names(trans_dt))) {
  stop("trans file is missing required columns: snps, gene")
}

snp_pattern <- "^(Chr[^_]+)_([0-9]+)$"
message("Parsing SNP coordinates from trans file (pattern: ", snp_pattern, ")...")
snp_dt <- unique(trans_dt[, .(snps)])
snp_dt <- snp_dt[grepl(snp_pattern, snps)]
snp_dt[, chr := sub(snp_pattern, "\\1", snps)]
snp_dt[, pos := as.integer(sub(snp_pattern, "\\2", snps))]
snp_dt <- snp_dt[!is.na(pos)]

if (nrow(snp_dt) == 0L) {
  stop("No valid SNP positions parsed from trans file.")
}

# Define 20 kb windows (±10 kb) around each SNP
snp_dt[, start := pmax(pos - 10000L, 1L)]
snp_dt[, end := pos + 10000L]
setnames(snp_dt, "chr", "seqid")

message("Loading gene annotations from GFF3 (gene features only)...")
gff_dt <- fread(
  gff_file,
  sep = "\t",
  header = FALSE,
  quote = "",
  col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
  data.table = TRUE
)

gene_dt <- gff_dt[type == "gene", .(seqid, start, end, attributes)]

if (nrow(gene_dt) == 0L) {
  stop("No 'gene' features found in GFF3 file.")
}

gene_dt[, gene_id := sub(".*ID=([^;]+);?.*", "\\1", attributes)]
gene_dt <- gene_dt[!is.na(gene_id)]

if (nrow(gene_dt) == 0L) {
  stop("No gene IDs could be parsed from GFF3 attributes.")
}

gene_iv <- gene_dt[, .(seqid, start, end, gene_id)]

message("Finding genes within ±10 kb of each SNP (20 kb window)...")
setkey(gene_iv, seqid, start, end)
setkey(snp_dt, seqid, start, end)

ov <- foverlaps(snp_dt, gene_iv, nomatch = 0L)

if (nrow(ov) == 0L) {
  stop("No genes found within 20 kb windows around the SNPs.")
}

cis_window <- unique(ov[, .(snps, window_gene = gene_id)])

message("Number of SNP–window-gene (20 kb) pairs: ", nrow(cis_window))

message("Building SNP–window-gene–trans triplets...")
setkey(cis_window, snps)
setkey(trans_dt, snps)

triplets <- cis_window[trans_dt, allow.cartesian = TRUE, nomatch = 0L]

if (nrow(triplets) == 0L) {
  stop("Join produced zero triplets; check that SNP IDs match between cis-window and trans tables.")
}

setnames(triplets, "gene", "trans_gene")

triplets_out <- triplets[
  ,
  .(snps, window_gene, trans_gene)
]

message("Constructed ", nrow(triplets_out), " SNP–cis-window–trans triplets (20 kb, winter). Writing to: ", out_file)
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
fwrite(triplets_out, out_file, sep = "\t")

message("Done (winter 20 kb triplets).")

