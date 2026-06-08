#!/usr/bin/env Rscript
# =============================================================================
# run_eqtl_pipeline.R
#
# Single-script eQTL pipeline:
#   Step 1. Align expression matrix, genotype matrix and covariate matrix on
#           the set of samples that are common to all three input files and
#           write the aligned matrices in Matrix eQTL format.
#   Step 2. Run Matrix eQTL (cis + trans) on the aligned matrices using the
#           supplied SNP location and gene position files.
#
# Inspired by:
#   ref_scripts/run7.R                  -> driver / argument layout
#   ref_scripts/cis_trans_matrix_eqtl_main.R -> Matrix eQTL invocation
#   ref_scripts/eqtl_config.R           -> parameter defaults
#
# Defaults are set so that running the script with no arguments analyses the
# `summer` dataset located at
#   050_eqtl/summer/data
# and writes everything to
#   050_eqtl/summer/results
# =============================================================================

suppressPackageStartupMessages({
  library(MatrixEQTL)
  library(data.table)
})

# -----------------------------------------------------------------------------
# Argument parsing (same `--key value` style as the reference scripts)
# -----------------------------------------------------------------------------
parse_args <- function(args) {
  out <- list()
  if (length(args) == 0) return(out)
  if (length(args) %% 2 != 0) {
    stop("Arguments must be supplied as `--key value` pairs.")
  }
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    out[[key]] <- args[i + 1]
  }
  out
}

params <- parse_args(commandArgs(trailingOnly = TRUE))

get_param <- function(name, default) {
  v <- params[[name]]
  if (is.null(v) || !nzchar(v)) default else v
}

# -----------------------------------------------------------------------------
# Default paths -- tuned for the `summer` dataset
# -----------------------------------------------------------------------------
base_dir    <- "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/050_eqtl/summer"
data_dir    <- get_param("data_dir",    file.path(base_dir, "data"))
results_dir <- get_param("results_dir", file.path(base_dir, "results"))
label       <- get_param("label",       "summer")

genotype_file     <- get_param("genotype_file",     file.path(data_dir, "genotype_matrix_output1.txt"))
expression_file   <- get_param("expression_file",   file.path(data_dir, "gene_count_matrix_summer_vsd_renamed.csv"))
covariate_file    <- get_param("covariate_file",    file.path(data_dir, "covariates_for_matrix_eqtl.txt"))
snp_location_file <- get_param("snp_location_file", file.path(data_dir, "snp_location_test1.txt"))
genepos_file      <- get_param("genepos_file",      file.path(data_dir, "genepos_clean.txt"))

pvOutputThreshold_cis   <- as.numeric(get_param("pv_cis",   "1e-5"))
pvOutputThreshold_trans <- as.numeric(get_param("pv_trans", "1e-5"))
cisDist                 <- as.numeric(get_param("cis_dist", "1e6"))
use_covariates          <- tolower(get_param("use_covariates", "TRUE")) %in% c("true", "t", "yes", "1")
skip_alignment          <- tolower(get_param("skip_alignment", "FALSE")) %in% c("true", "t", "yes", "1")
expression_sep          <- get_param("expression_sep", ",")  # csv input by default
genotype_sep            <- get_param("genotype_sep",   "\t")
covariate_sep           <- get_param("covariate_sep",  "\t")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Derived (aligned) file paths
geno_aligned <- file.path(results_dir, paste0("genotype_matrix_aligned_",   label, ".txt"))
expr_aligned <- file.path(results_dir, paste0("expression_matrix_aligned_", label, ".txt"))
cov_aligned  <- file.path(results_dir, paste0("covariates_matrix_aligned_", label, ".txt"))

# Matrix eQTL output paths
cis_out_file   <- file.path(results_dir, paste0("cis_eqtl_results_",   label, ".txt"))
trans_out_file <- file.path(results_dir, paste0("trans_eqtl_results_", label, ".txt"))
log_file       <- file.path(results_dir, paste0("matrix_eqtl_log_",    label, ".txt"))
hist_file      <- file.path(results_dir, paste0("matrix_eqtl_pvalue_histogram_", label, ".png"))
samples_file   <- file.path(results_dir, paste0("common_samples_",     label, ".txt"))

log_msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")

log_msg("==================== eQTL pipeline started ====================")
log_msg("Data dir         : ", data_dir)
log_msg("Results dir      : ", results_dir)
log_msg("Label            : ", label)
log_msg("Use covariates   : ", use_covariates)
log_msg("Skip alignment   : ", skip_alignment)
log_msg("pv cis / trans   : ", pvOutputThreshold_cis, " / ", pvOutputThreshold_trans)
log_msg("cis distance (bp): ", cisDist)

# =============================================================================
# Step 1: ALIGN expression, genotype and covariate matrices on common samples
# =============================================================================
if (skip_alignment &&
    file.exists(geno_aligned) &&
    file.exists(expr_aligned) &&
    (!use_covariates || file.exists(cov_aligned))) {
  log_msg("skip_alignment=TRUE and aligned files already exist -- skipping Step 1.")
} else {
  log_msg("---- Step 1: aligning samples ----")

  # --- Read headers only, to determine sample sets ---------------------------
  geno_samples <- unlist(strsplit(readLines(genotype_file,   n = 1), genotype_sep,   fixed = TRUE))[-1]
  expr_samples <- unlist(strsplit(readLines(expression_file, n = 1), expression_sep, fixed = TRUE))[-1]

  # Covariate file: header has NO row-name label, so every header field is a
  # sample name (do NOT drop the first field here, unlike genotype/expression).
  cov_header_line <- readLines(covariate_file, n = 1)
  cov_samples     <- unlist(strsplit(cov_header_line, covariate_sep, fixed = TRUE))

  log_msg("Samples - genotype : ", length(geno_samples))
  log_msg("Samples - expression: ", length(expr_samples))
  log_msg("Samples - covariates: ", length(cov_samples))

  sample_lists <- list(geno_samples, expr_samples)
  if (use_covariates) sample_lists <- c(sample_lists, list(cov_samples))
  common_samples <- Reduce(intersect, sample_lists)
  # Preserve order from the genotype file (most expensive to re-sort later)
  common_samples <- geno_samples[geno_samples %in% common_samples]

  if (length(common_samples) == 0) {
    stop("No samples are common to all input files. Check sample naming.")
  }
  log_msg("Common samples   : ", length(common_samples))
  writeLines(common_samples, samples_file)

  # --- Genotype matrix (the big one: ~37 GB) ---------------------------------
  log_msg("Loading aligned genotype matrix via data.table::fread (this may take a while)...")
  geno_first_col <- unlist(strsplit(readLines(genotype_file, n = 1), genotype_sep, fixed = TRUE))[1]
  geno_dt <- fread(
    genotype_file,
    sep    = genotype_sep,
    header = TRUE,
    select = c(geno_first_col, common_samples),
    check.names = FALSE,
    showProgress = TRUE
  )
  setcolorder(geno_dt, c(geno_first_col, common_samples))
  log_msg("Writing aligned genotype matrix: ", geno_aligned)
  fwrite(geno_dt, file = geno_aligned, sep = "\t", quote = FALSE, na = "NA")
  rm(geno_dt); invisible(gc())

  # --- Expression matrix -----------------------------------------------------
  log_msg("Loading aligned expression matrix...")
  expr_first_col <- unlist(strsplit(readLines(expression_file, n = 1), expression_sep, fixed = TRUE))[1]
  expr_dt <- fread(
    expression_file,
    sep    = expression_sep,
    header = TRUE,
    select = c(expr_first_col, common_samples),
    check.names = FALSE,
    showProgress = TRUE
  )
  setcolorder(expr_dt, c(expr_first_col, common_samples))
  log_msg("Writing aligned expression matrix: ", expr_aligned)
  fwrite(expr_dt, file = expr_aligned, sep = "\t", quote = FALSE, na = "NA")
  rm(expr_dt); invisible(gc())

  # --- Covariate matrix ------------------------------------------------------
  if (use_covariates) {
    log_msg("Loading and aligning covariate matrix...")
    cov_mat <- read.table(
      covariate_file,
      header = TRUE,
      sep    = covariate_sep,
      row.names = 1,         # 1st column of every data row is the covariate name
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    missing_cov_samples <- setdiff(common_samples, colnames(cov_mat))
    if (length(missing_cov_samples) > 0) {
      stop("Covariate matrix is missing samples that survived intersection: ",
           paste(head(missing_cov_samples, 5), collapse = ", "))
    }
    cov_mat <- cov_mat[, common_samples, drop = FALSE]
    cov_out <- data.table(id = rownames(cov_mat), cov_mat, check.names = FALSE)
    log_msg("Writing aligned covariate matrix: ", cov_aligned)
    fwrite(cov_out, file = cov_aligned, sep = "\t", quote = FALSE, na = "NA")
    rm(cov_mat, cov_out); invisible(gc())
  } else {
    log_msg("Covariates disabled -- not writing covariate file.")
  }
}

# =============================================================================
# Step 2: RUN Matrix eQTL on the aligned matrices
# =============================================================================
log_msg("---- Step 2: running Matrix eQTL ----")

useModel <- modelLINEAR

# --- Genotype (SNP) SlicedData -----------------------------------------------
snps <- SlicedData$new()
snps$fileDelimiter      <- "\t"
snps$fileOmitCharacters <- "NA"
snps$fileSkipRows       <- 1
snps$fileSkipColumns    <- 1
snps$fileSliceSize      <- 2000
log_msg("Loading aligned genotypes into SlicedData: ", geno_aligned)
snps$LoadFile(geno_aligned)

# --- Expression SlicedData ---------------------------------------------------
gene <- SlicedData$new()
gene$fileDelimiter      <- "\t"
gene$fileOmitCharacters <- "NA"
gene$fileSkipRows       <- 1
gene$fileSkipColumns    <- 1
gene$fileSliceSize      <- 2000
log_msg("Loading aligned expression into SlicedData: ", expr_aligned)
gene$LoadFile(expr_aligned)

# --- Covariate SlicedData ----------------------------------------------------
cvrt <- SlicedData$new()
if (use_covariates) {
  cvrt$fileDelimiter      <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows       <- 1
  cvrt$fileSkipColumns    <- 1
  cvrt$fileSliceSize      <- 2000
  log_msg("Loading aligned covariates into SlicedData: ", cov_aligned)
  cvrt$LoadFile(cov_aligned)
}

# --- SNP and gene position tables --------------------------------------------
log_msg("Loading SNP location file: ", snp_location_file)
snpspos <- read.table(snp_location_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
# Expected columns: SNP, CHR, POS (in that order)
required_snp_cols <- c("SNP", "CHR", "POS")
if (!all(required_snp_cols %in% colnames(snpspos))) {
  stop("SNP location file must contain columns: ", paste(required_snp_cols, collapse = ", "))
}
snpspos <- snpspos[, required_snp_cols]
snpspos$CHR <- as.character(snpspos$CHR)
snpspos$POS <- as.numeric(snpspos$POS)

log_msg("Loading gene position file: ", genepos_file)
genepos <- read.table(genepos_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
required_gene_cols <- c("geneid", "chr", "left", "right")
if (!all(required_gene_cols %in% colnames(genepos))) {
  stop("Gene position file must contain columns: ", paste(required_gene_cols, collapse = ", "))
}
genepos <- genepos[, required_gene_cols]
genepos$chr   <- as.character(genepos$chr)
genepos$left  <- as.numeric(genepos$left)
genepos$right <- as.numeric(genepos$right)

# --- Matrix eQTL -------------------------------------------------------------
log_msg("Running Matrix_eQTL_main ...")
me <- Matrix_eQTL_main(
  snps                  = snps,
  gene                  = gene,
  cvrt                  = cvrt,
  output_file_name      = trans_out_file,
  pvOutputThreshold     = pvOutputThreshold_trans,
  useModel              = useModel,
  errorCovariance       = numeric(),
  verbose               = TRUE,
  output_file_name.cis  = cis_out_file,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos               = snpspos,
  genepos               = genepos,
  cisDist               = cisDist,
  pvalue.hist           = TRUE,
  min.pv.by.genesnp     = FALSE,
  noFDRsaveMemory       = FALSE
)

# Re-write the eqtl tables explicitly to guarantee full FDR columns are saved.
write.table(me$cis$eqtls,   file = cis_out_file,   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(me$trans$eqtls, file = trans_out_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Log + histogram
sink(log_file)
cat("Analysis completed in: ", me$time.in.sec, " seconds\n", sep = "")
cat("cis  eQTLs tested: ", me$cis$ntests,  "  detected (<= threshold): ", nrow(me$cis$eqtls),  "\n", sep = "")
cat("trans eQTLs tested: ", me$trans$ntests, "  detected (<= threshold): ", nrow(me$trans$eqtls), "\n", sep = "")
sink()

png(hist_file, width = 800, height = 600)
plot(me)
dev.off()

log_msg("Matrix eQTL finished in ", round(me$time.in.sec, 1), " seconds.")
log_msg("cis  results : ", cis_out_file)
log_msg("trans results: ", trans_out_file)
log_msg("log file     : ", log_file)
log_msg("histogram    : ", hist_file)
log_msg("==================== eQTL pipeline done ====================")
