#!/usr/bin/env Rscript
#############################################################################
####### Extract Winter-Biased Genes from DESeq2 Results (No Threshold) #####
#############################################################################
# This script extracts winter-biased genes from DESeq2 results (no filtering version)
# Criteria: FDR < 0.05 and log2FC > 0 (winter-biased = higher in winter)

setwd("/mnt/winterprojectceph/new_winter/001_deseq/test")

# Output directory
output_dir <- "results/without_threshold"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#===============================================================================
# Load DESeq2 results
#===============================================================================
cat("Extracting winter-biased genes (no threshold version)...\n")

# Load the DESeq2 results file (FDR < 0.001, no threshold version)
deseq2_results_file <- file.path(output_dir, "DESeq2results_winter_VS_summer_filtered_FDR0.001_no_threshold.txt")

if (!file.exists(deseq2_results_file)) {
  stop("DESeq2 results file not found: ", deseq2_results_file, "\nPlease run deseq2_summer_vs_winter_no_threshold.R first.")
}

# Read DESeq2 results
deseq2_results <- read.table(deseq2_results_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Loaded DESeq2 results:", nrow(deseq2_results), "genes\n")

#===============================================================================
# Extract winter-biased genes
#===============================================================================
# Winter-biased genes: FDR < 0.05 and log2FC > 0 (higher expression in winter)
alpha <- 0.05

winter_biased <- deseq2_results[
  !is.na(deseq2_results$padj) & 
  deseq2_results$padj < alpha & 
  deseq2_results$log2FoldChange > 0,
]

cat("Winter-biased genes (FDR < ", alpha, ", log2FC > 0): ", nrow(winter_biased), "\n", sep="")

# Rename gene_id column if needed (some files use different column names)
if (!"gene_id" %in% colnames(winter_biased)) {
  if ("gene" %in% colnames(winter_biased)) {
    colnames(winter_biased)[colnames(winter_biased) == "gene"] <- "gene_id"
  } else if (nrow(winter_biased) > 0 && "X" %in% colnames(winter_biased)) {
    # If gene IDs are in row names, add as column
    winter_biased$gene_id <- rownames(winter_biased)
  }
}

# Save winter-biased genes (FDR < 0.05, log2FC > 0)
output_file <- file.path(output_dir, "winter_biased_genes_FDR0.05_no_threshold.txt")
write.table(winter_biased, file = output_file, 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Saved winter-biased genes file: ", output_file, "\n")

#===============================================================================
# Extract strict winter-biased genes (log2FC >= 1, i.e., at least 2-fold higher)
#===============================================================================
winter_biased_strict <- winter_biased[winter_biased$log2FoldChange >= 1, ]
cat("Winter-biased genes (log2FC >= 1): ", nrow(winter_biased_strict), "\n", sep="")

# Save strict winter-biased genes
output_file_strict <- file.path(output_dir, "winter_biased_genes_FDR0.05_log2FC1_no_threshold.txt")
write.table(winter_biased_strict, file = output_file_strict, 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Saved strict winter-biased genes file (log2FC >= 1): ", output_file_strict, "\n")

cat("\n===============================================================================\n")
cat("EXTRACTION COMPLETE\n")
cat("===============================================================================\n")
cat("Winter-biased genes (FDR < 0.05, log2FC > 0): ", nrow(winter_biased), "\n")
cat("Winter-biased genes (FDR < 0.05, log2FC >= 1): ", nrow(winter_biased_strict), "\n")
cat("===============================================================================\n")
