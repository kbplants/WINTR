#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# Usage:
#   Rscript scripts/run_mediation_20kb_winter.R .
# or:
#   Rscript scripts/run_mediation_20kb_winter.R /path/to/032_cis_mediation_fdr20
#
# Inputs:
#   - results/winter/cis_trans_triplets_all_20kb_winter.tsv  (SNP–window-gene–trans triplets)
#   - data/wint_raw_cnt_mtx/gene_count_matrix_winter.csv  (winter expression matrix)
#   - data/genotype_matrix_output1.txt
# Output:
#   - results/winter/mediation_results_20kb.tsv
#
# This script reuses the same mediation regression logic as run_mediation.R,
# but operates on pre-built 20 kb window triplets instead of cis/trans eQTL tables.

args <- commandArgs(trailingOnly = TRUE)

base_dir <- if (length(args) >= 1) args[[1]] else "."

triplet_file <- file.path(base_dir, "results", "winter", "cis_trans_triplets_all_20kb_winter.tsv")
expr_file <- file.path(base_dir, "data", "wint_raw_cnt_mtx", "gene_count_matrix_winter.csv")
geno_file <- file.path(base_dir, "data", "genotype_matrix_output1.txt")
out_file <- file.path(base_dir, "results", "winter", "mediation_results_20kb.tsv")

message("Loading 20 kb window SNP–window-gene–trans triplets (winter)...")
triplets_dt <- fread(triplet_file)

required_cols <- c("snps", "window_gene", "trans_gene")
if (!all(required_cols %in% names(triplets_dt))) {
  stop("Triplet file is missing required columns: ", paste(setdiff(required_cols, names(triplets_dt)), collapse = ", "))
}

overlap_snps <- unique(triplets_dt$snps)
if (length(overlap_snps) == 0L) {
  stop("No SNPs found in 20 kb triplets file.")
}

message("Loading expression matrix (winter; this may take a while)...")
expr_dt <- fread(expr_file)
if (!"gene_id" %in% names(expr_dt)) {
  stop("Expression matrix must have a 'gene_id' column.")
}
gene_ids <- expr_dt$gene_id
expr_mat <- as.matrix(expr_dt[, -"gene_id", with = FALSE])
rownames(expr_mat) <- gene_ids

sample_ids_expr <- colnames(expr_mat)

message("Loading genotype matrix (this may take a while)...")
geno_dt <- fread(geno_file)

first_col_name <- names(geno_dt)[1L]
if (!first_col_name %in% c("snps", "SNP", "snp_id")) {
  stop("First column of genotype matrix must be SNP identifier (snps/SNP/snp_id).")
}

setnames(geno_dt, first_col_name, "snps")

geno_snps <- geno_dt$snps
geno_mat <- as.matrix(geno_dt[, -"snps", with = FALSE])
sample_ids_geno <- colnames(geno_mat)

common_samples <- intersect(sample_ids_expr, sample_ids_geno)
if (length(common_samples) < 10L) {
  stop("Too few overlapping samples between genotype and expression matrices.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
geno_mat <- geno_mat[, common_samples, drop = FALSE]

message("Preparing SNP index lookup (winter 20 kb)...")
geno_index <- match(overlap_snps, geno_snps)
names(geno_index) <- overlap_snps

if (any(is.na(geno_index))) {
  missing_snps <- names(geno_index)[is.na(geno_index)]
  warning(length(missing_snps), " SNPs from triplets missing in genotype matrix; they will be skipped.")
}

message("Starting mediation regressions for ", nrow(triplets_dt), " winter 20 kb triplets...")

run_mediation_for_triplet <- function(i) {
  snp_id <- triplets_dt$snps[i]
  window_gene <- triplets_dt$window_gene[i]
  trans_gene <- triplets_dt$trans_gene[i]

  gi <- geno_index[[snp_id]]
  if (is.na(gi)) return(NULL)

  if (!window_gene %in% rownames(expr_mat) || !trans_gene %in% rownames(expr_mat)) {
    return(NULL)
  }

  x_snp <- as.numeric(geno_mat[gi, ])
  x_cis <- as.numeric(expr_mat[window_gene, ])
  y_trans <- as.numeric(expr_mat[trans_gene, ])

  if (sd(x_snp, na.rm = TRUE) == 0 || sd(x_cis, na.rm = TRUE) == 0 || sd(y_trans, na.rm = TRUE) == 0) {
    return(NULL)
  }

  complete_idx <- which(!is.na(x_snp) & !is.na(x_cis) & !is.na(y_trans))
  if (length(complete_idx) < 10L) return(NULL)

  x_snp <- x_snp[complete_idx]
  x_cis <- x_cis[complete_idx]
  y_trans <- y_trans[complete_idx]

  cis_fit <- try(lm(x_cis ~ x_snp), silent = TRUE)
  if (inherits(cis_fit, "try-error")) return(NULL)
  cis_coef <- summary(cis_fit)$coef
  a <- cis_coef["x_snp", "Estimate"]
  se_a <- cis_coef["x_snp", "Std. Error"]

  total_fit <- try(lm(y_trans ~ x_snp), silent = TRUE)
  if (inherits(total_fit, "try-error")) return(NULL)
  total_coef <- summary(total_fit)$coef
  b_unadj <- total_coef["x_snp", "Estimate"]
  se_b_unadj <- total_coef["x_snp", "Std. Error"]

  med_fit <- try(lm(y_trans ~ x_snp + x_cis), silent = TRUE)
  if (inherits(med_fit, "try-error")) return(NULL)
  med_coef <- summary(med_fit)$coef
  if (!"x_snp" %in% rownames(med_coef) || !"x_cis" %in% rownames(med_coef)) return(NULL)
  b_adj <- med_coef["x_snp", "Estimate"]
  se_b_adj <- med_coef["x_snp", "Std. Error"]
  b <- med_coef["x_cis", "Estimate"]
  se_b <- med_coef["x_cis", "Std. Error"]

  if (abs(b_unadj) < 1e-8) {
    mp <- NA_real_
  } else {
    mp <- (abs(b_unadj) - abs(b_adj)) / abs(b_unadj)
  }

  ab <- a * b
  se_ab <- sqrt((b^2) * (se_a^2) + (a^2) * (se_b^2))
  z_sobel <- if (se_ab > 0) ab / se_ab else NA_real_
  p_sobel <- if (!is.na(z_sobel)) 2 * pnorm(-abs(z_sobel)) else NA_real_

  data.table(
    snps = snp_id,
    window_gene = window_gene,
    trans_gene = trans_gene,
    a = a,
    se_a = se_a,
    b_unadj = b_unadj,
    se_b_unadj = se_b_unadj,
    b_adj = b_adj,
    se_b_adj = se_b_adj,
    b_mediator = b,
    se_b_mediator = se_b,
    mediation_prop = mp,
    indirect_effect = ab,
    se_indirect = se_ab,
    z_sobel = z_sobel,
    p_sobel = p_sobel,
    n_samples = length(complete_idx)
  )
}

results_list <- vector("list", nrow(triplets_dt))
for (i in seq_len(nrow(triplets_dt))) {
  if (i %% 1000L == 0L) {
    message("Processed ", i, " / ", nrow(triplets_dt), " winter 20 kb triplets...")
  }
  res <- run_mediation_for_triplet(i)
  results_list[[i]] <- res
}

mediation_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
mediation_dt <- mediation_dt[!is.na(a) & !is.na(b_unadj) & !is.na(b_adj)]

message("Writing winter 20 kb mediation results to: ", out_file)
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
fwrite(mediation_dt, out_file, sep = "\t")

message("Done (winter 20 kb mediation).")

