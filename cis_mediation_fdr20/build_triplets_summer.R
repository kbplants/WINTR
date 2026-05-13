#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# Usage:
#   Rscript scripts/build_triplets_summer.R .
# or:
#   Rscript scripts/build_triplets_summer.R /path/to/032_cis_mediation_fdr20
#
# Inputs (in data/): cis_eqtl_summer_xylem_gene_filter_pf_all_new_fdr_20.txt,
#                    trans_eqtl_results_summer_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt
# Output: results/summer/cis_trans_triplets_all.tsv

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "."

cis_file <- file.path(base_dir, "data", "cis_eqtl_summer_xylem_gene_filter_pf_all_new_fdr_20.txt")
trans_file <- file.path(base_dir, "data", "trans_eqtl_results_summer_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt")
out_file <- file.path(base_dir, "results", "summer", "cis_trans_triplets_all.tsv")

message("Loading cis and trans eQTL results (summer)...")
cis_dt <- fread(cis_file)
trans_dt <- fread(trans_file)

required_cols <- c("snps", "gene", "statistic", "pvalue", "FDR", "beta")
if (!all(required_cols %in% names(cis_dt))) {
  stop("cis file is missing required columns: ", paste(setdiff(required_cols, names(cis_dt)), collapse = ", "))
}
if (!all(required_cols %in% names(trans_dt))) {
  stop("trans file is missing required columns: ", paste(setdiff(required_cols, names(trans_dt)), collapse = ", "))
}

message("Determining overlapping SNPs for cis and trans (summer, no additional FDR filter)...")
overlap_snps <- intersect(unique(cis_dt$snps), unique(trans_dt$snps))
cis_dt <- cis_dt[snps %in% overlap_snps]
trans_dt <- trans_dt[snps %in% overlap_snps]

if (nrow(cis_dt) == 0L || nrow(trans_dt) == 0L) {
  stop("No overlapping SNPs between cis and trans tables after restricting to shared SNPs.")
}

message("Building all cis–trans gene triplets per overlapping SNP (summer)...")
setkey(cis_dt, snps)
setkey(trans_dt, snps)

# Cartesian join over genes within each SNP (no p/FDR thresholding)
triplets <- cis_dt[trans_dt, allow.cartesian = TRUE, nomatch = 0L]

if (nrow(triplets) == 0L) {
  stop("Join produced zero triplets; check that cis and trans files share SNP IDs.")
}

# Rename and select columns for clarity
setnames(triplets,
         old = c("gene", "statistic", "pvalue", "FDR", "beta",
                 "i.gene", "i.statistic", "i.pvalue", "i.FDR", "i.beta"),
         new = c("cis_gene", "cis_statistic", "cis_pvalue", "cis_FDR", "cis_beta",
                 "trans_gene", "trans_statistic", "trans_pvalue", "trans_FDR", "trans_beta"))

triplets_out <- triplets[
  ,
  .(snps,
    cis_gene, cis_beta, cis_statistic, cis_pvalue, cis_FDR,
    trans_gene, trans_beta, trans_statistic, trans_pvalue, trans_FDR)
]

message("Constructed ", nrow(triplets_out), " SNP–cis–trans triplets (summer). Writing to: ", out_file)
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
fwrite(triplets_out, out_file, sep = "\t")

message("Done (summer triplets).")

