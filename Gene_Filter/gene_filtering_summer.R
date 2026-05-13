#!/usr/bin/env Rscript
# =============================================================================
# Script: gene_filtering_summer.R
# Purpose: Gene filtering pipeline for summer gene count matrix
# Author: Generated for winter project
# Date: 2025-01-02
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(readxl)
})

# =============================================================================
# Configuration
# =============================================================================

# Base directory
base_dir <- "/mnt/winterprojectceph/winter_pipe/combine_batch_results"
scripts_dir <- file.path(base_dir, "gf_cnt_mtx", "scripts")
results_dir <- file.path(base_dir, "gf_cnt_mtx", "results")
raw_cnt_dir <- file.path(base_dir, "raw_cnt_mtx", "results")

# Input file
input_file <- file.path(raw_cnt_dir, "gene_count_matrix_summer.csv")

# Gene stats file (optional - for min expression filtering)
gene_stats_file <- "~/Downloads/gene_stats_wint.xlsx"

# Output files
output_filtered <- file.path(results_dir, "gene_count_matrix_summer_filtered.csv")
output_vsd <- file.path(results_dir, "gene_count_matrix_summer_vsd.csv")
output_summary <- file.path(results_dir, "gene_filtering_summary.txt")

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Open summary file for writing
summary_file <- file(output_summary, "w")

cat(paste(rep("=", 60), collapse = ""), "\n", file = summary_file)
cat("Gene Filtering Pipeline for Summer Count Matrix\n", file = summary_file)
cat(paste(rep("=", 60), collapse = ""), "\n\n", file = summary_file)
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Gene Filtering Pipeline for Summer Count Matrix\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# =============================================================================
# Step 1: Load gene count matrix
# =============================================================================

cat("Step 1: Loading gene count matrix...\n")
cat("Step 1: Loading gene count matrix...\n", file = summary_file)

cat(sprintf("  Reading: %s\n", input_file))
cat(sprintf("  Reading: %s\n", input_file), file = summary_file)

# Read gene count matrix
mtx <- read.csv(input_file, row.names = 1, check.names = FALSE)

cat(sprintf("  Original matrix: %d genes × %d samples\n", nrow(mtx), ncol(mtx)))
cat(sprintf("  Original matrix: %d genes × %d samples\n", nrow(mtx), ncol(mtx)), file = summary_file)

# =============================================================================
# Step 2: Initial filtering - remove zeros
# =============================================================================

cat("\nStep 2: Initial filtering - removing zero rows and columns...\n")
cat("\nStep 2: Initial filtering - removing zero rows and columns...\n", file = summary_file)

# Replace NAs with 0
mtx[is.na(mtx)] <- 0

# Remove all rows where sum values are 0
mtx_filtered <- mtx[rowSums(mtx) != 0, ]

# Remove all columns where sum values are 0
mtx_filtered <- mtx_filtered[, colSums(mtx_filtered) != 0]

cat(sprintf("  After zero removal: %d genes × %d samples\n", nrow(mtx_filtered), ncol(mtx_filtered)))
cat(sprintf("  After zero removal: %d genes × %d samples\n", nrow(mtx_filtered), ncol(mtx_filtered)), file = summary_file)
cat(sprintf("  Removed: %d genes, %d samples\n", 
            nrow(mtx) - nrow(mtx_filtered), 
            ncol(mtx) - ncol(mtx_filtered)))
cat(sprintf("  Removed: %d genes, %d samples\n", 
            nrow(mtx) - nrow(mtx_filtered), 
            ncol(mtx) - ncol(mtx_filtered)), file = summary_file)

# =============================================================================
# Step 3: IQR-based sample filtering
# =============================================================================

cat("\nStep 3: IQR-based sample filtering...\n")
cat("\nStep 3: IQR-based sample filtering...\n", file = summary_file)

# Calculate number of genes expressed in every sample
genes_expressed_per_sample <- colSums(mtx_filtered > 0)

# Calculate IQR-based threshold
Q1 <- quantile(genes_expressed_per_sample, 0.25)
Q3 <- quantile(genes_expressed_per_sample, 0.75)
IQR_val <- IQR(genes_expressed_per_sample)

# Set threshold (can be adjusted)
threshold <- Q1 - 1.5 * IQR_val

# Print values
cat(sprintf("  Q1: %.2f, Q3: %.2f, IQR: %.2f, Threshold: %.2f\n", Q1, Q3, IQR_val, threshold))
cat(sprintf("  Q1: %.2f, Q3: %.2f, IQR: %.2f, Threshold: %.2f\n", Q1, Q3, IQR_val, threshold), file = summary_file)

# Retain samples which pass the threshold
samples_to_keep <- names(genes_expressed_per_sample[genes_expressed_per_sample > threshold])
samples_removed <- setdiff(colnames(mtx_filtered), samples_to_keep)

cat(sprintf("  Samples to keep: %d\n", length(samples_to_keep)))
cat(sprintf("  Samples removed: %d\n", length(samples_removed)))
cat(sprintf("  Samples to keep: %d\n", length(samples_to_keep)), file = summary_file)
cat(sprintf("  Samples removed: %d\n", length(samples_removed)), file = summary_file)

# =============================================================================
# Step 4: IQR-based gene filtering
# =============================================================================

cat("\nStep 4: IQR-based gene filtering...\n")
cat("\nStep 4: IQR-based gene filtering...\n", file = summary_file)

# Count how many samples each gene is expressed in
samples_expressed_per_gene <- rowSums(mtx_filtered > 0)

# Calculate IQR-based threshold
Q1_gene <- quantile(samples_expressed_per_gene, 0.25)
Q3_gene <- quantile(samples_expressed_per_gene, 0.75)
IQR_gene <- IQR(samples_expressed_per_gene)

# Set threshold (adjustable)
threshold_gene <- Q1_gene - 1.5 * IQR_gene

# Print values for sanity check
cat(sprintf("  Q1: %.2f, Q3: %.2f, IQR: %.2f, Threshold: %.2f\n", 
            Q1_gene, Q3_gene, IQR_gene, threshold_gene))
cat(sprintf("  Q1: %.2f, Q3: %.2f, IQR: %.2f, Threshold: %.2f\n", 
            Q1_gene, Q3_gene, IQR_gene, threshold_gene), file = summary_file)

# Keep genes expressed in enough samples
genes_to_keep <- names(samples_expressed_per_gene[samples_expressed_per_gene > threshold_gene])
genes_removed <- setdiff(rownames(mtx_filtered), genes_to_keep)

cat(sprintf("  Genes to keep: %d\n", length(genes_to_keep)))
cat(sprintf("  Genes removed: %d\n", length(genes_removed)))
cat(sprintf("  Genes to keep: %d\n", length(genes_to_keep)), file = summary_file)
cat(sprintf("  Genes removed: %d\n", length(genes_removed)), file = summary_file)

# Apply both filters
filtered_mtx <- mtx_filtered[genes_to_keep, samples_to_keep]

cat(sprintf("  After IQR filtering: %d genes × %d samples\n", 
            nrow(filtered_mtx), ncol(filtered_mtx)))
cat(sprintf("  After IQR filtering: %d genes × %d samples\n", 
            nrow(filtered_mtx), ncol(filtered_mtx)), file = summary_file)

# =============================================================================
# Step 5: Filter by zero percentage (0-10%)
# =============================================================================

cat("\nStep 5: Filtering genes by zero percentage (0-10%)...\n")
cat("\nStep 5: Filtering genes by zero percentage (0-10%)...\n", file = summary_file)

# Filter according to 0-10% zero percentage
zero_percent <- rowMeans(filtered_mtx == 0) * 100

expr_data <- data.frame(
  gene_id = rownames(filtered_mtx),
  zero_percent = zero_percent
)

expr_data <- expr_data %>%
  mutate(zero_bin = cut(zero_percent,
                        breaks = seq(0, 100, by = 10),
                        include.lowest = TRUE,
                        right = FALSE,
                        labels = paste(seq(0, 90, by = 10), seq(10, 100, by = 10), sep = "-")))

expr_data_filtered <- expr_data %>% filter(zero_bin == "0-10")
filtered_mtx_final <- filtered_mtx[expr_data_filtered$gene_id, ]

cat(sprintf("  Genes with 0-10%% zeros: %d\n", nrow(filtered_mtx_final)))
cat(sprintf("  Genes removed: %d\n", nrow(filtered_mtx) - nrow(filtered_mtx_final)))
cat(sprintf("  Genes with 0-10%% zeros: %d\n", nrow(filtered_mtx_final)), file = summary_file)
cat(sprintf("  Genes removed: %d\n", nrow(filtered_mtx) - nrow(filtered_mtx_final)), file = summary_file)

# =============================================================================
# Step 6: Filter by 10th percentile of average expression
# =============================================================================

cat("\nStep 6: Filtering by 10th percentile of average expression...\n")
cat("\nStep 6: Filtering by 10th percentile of average expression...\n", file = summary_file)

# Calculate 10th percentile threshold
gene_avg_expr <- rowMeans(filtered_mtx_final, na.rm = TRUE)
threshold_10th <- quantile(gene_avg_expr, probs = 0.10, na.rm = TRUE)

cat(sprintf("  10th percentile threshold: %.4f\n", threshold_10th))
cat(sprintf("  10th percentile threshold: %.4f\n", threshold_10th), file = summary_file)

# Filter genes above 10th percentile
filtered_mtx1 <- filtered_mtx_final[gene_avg_expr > threshold_10th, ]

cat(sprintf("  Genes above 10th percentile: %d\n", nrow(filtered_mtx1)))
cat(sprintf("  Genes removed: %d\n", nrow(filtered_mtx_final) - nrow(filtered_mtx1)))
cat(sprintf("  Genes above 10th percentile: %d\n", nrow(filtered_mtx1)), file = summary_file)
cat(sprintf("  Genes removed: %d\n", nrow(filtered_mtx_final) - nrow(filtered_mtx1)), file = summary_file)

# =============================================================================
# Step 7: Filter by minimum expression > 15 (optional)
# =============================================================================

cat("\nStep 7: Filtering by minimum expression > 15...\n")
cat("\nStep 7: Filtering by minimum expression > 15...\n", file = summary_file)

# Check if gene stats file exists
if (file.exists(gene_stats_file)) {
  cat(sprintf("  Loading gene stats from: %s\n", gene_stats_file))
  cat(sprintf("  Loading gene stats from: %s\n", gene_stats_file), file = summary_file)
  
  expr_data_stats <- read_excel(gene_stats_file)
  min_gene_list <- expr_data_stats %>% 
    filter(min_expr > 15) %>% 
    pull(GeneID)
  
  # Filter to keep only genes in the list
  filtered_mtx_min <- filtered_mtx1[rownames(filtered_mtx1) %in% min_gene_list, ]
  
  cat(sprintf("  Genes with min expression > 15: %d\n", nrow(filtered_mtx_min)))
  cat(sprintf("  Genes removed: %d\n", nrow(filtered_mtx1) - nrow(filtered_mtx_min)))
  cat(sprintf("  Genes with min expression > 15: %d\n", nrow(filtered_mtx_min)), file = summary_file)
  cat(sprintf("  Genes removed: %d\n", nrow(filtered_mtx1) - nrow(filtered_mtx_min)), file = summary_file)
  
  # Use the min-filtered matrix for next steps
  filtered_mtx_final_step <- filtered_mtx_min
} else {
  cat(sprintf("  WARNING: Gene stats file not found: %s\n", gene_stats_file))
  cat(sprintf("  Skipping min expression filter, using previous filtered matrix\n"))
  cat(sprintf("  WARNING: Gene stats file not found: %s\n", gene_stats_file), file = summary_file)
  cat(sprintf("  Skipping min expression filter, using previous filtered matrix\n"), file = summary_file)
  
  filtered_mtx_final_step <- filtered_mtx1
}

cat(sprintf("  Final filtered matrix: %d genes × %d samples\n", 
            nrow(filtered_mtx_final_step), ncol(filtered_mtx_final_step)))
cat(sprintf("  Final filtered matrix: %d genes × %d samples\n", 
            nrow(filtered_mtx_final_step), ncol(filtered_mtx_final_step)), file = summary_file)

# =============================================================================
# Step 8: DESeq2 normalization and variance stabilization
# =============================================================================

cat("\nStep 8: DESeq2 normalization and variance stabilization...\n")
cat("\nStep 8: DESeq2 normalization and variance stabilization...\n", file = summary_file)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = filtered_mtx_final_step,
  colData = DataFrame(row.names = colnames(filtered_mtx_final_step)),
  design = ~ 1
)

# Estimate size factors
dds <- estimateSizeFactors(dds, type = "poscounts")

# Variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsd_matrix <- assay(vsd)

cat(sprintf("  VSD matrix: %d genes × %d samples\n", 
            nrow(vsd_matrix), ncol(vsd_matrix)))
cat(sprintf("  VSD matrix: %d genes × %d samples\n", 
            nrow(vsd_matrix), ncol(vsd_matrix)), file = summary_file)

# =============================================================================
# Step 9: Save results
# =============================================================================

cat("\nStep 9: Saving results...\n")
cat("\nStep 9: Saving results...\n", file = summary_file)

# Save filtered count matrix
filtered_mtx_out <- filtered_mtx_final_step
filtered_mtx_out$gene_id <- rownames(filtered_mtx_out)
filtered_mtx_out <- filtered_mtx_out[, c("gene_id", setdiff(colnames(filtered_mtx_out), "gene_id"))]

cat(sprintf("  Writing filtered matrix to: %s\n", output_filtered))
cat(sprintf("  Writing filtered matrix to: %s\n", output_filtered), file = summary_file)
write.csv(filtered_mtx_out, output_filtered, row.names = FALSE)
cat("  ✓ Filtered count matrix saved\n")
cat("  ✓ Filtered count matrix saved\n", file = summary_file)

# Save VSD matrix
vsd_mtx_out <- as.data.frame(vsd_matrix)
vsd_mtx_out$gene_id <- rownames(vsd_mtx_out)
vsd_mtx_out <- vsd_mtx_out[, c("gene_id", setdiff(colnames(vsd_mtx_out), "gene_id"))]

cat(sprintf("  Writing VSD matrix to: %s\n", output_vsd))
cat(sprintf("  Writing VSD matrix to: %s\n", output_vsd), file = summary_file)
write.csv(vsd_mtx_out, output_vsd, row.names = FALSE)
cat("  ✓ VSD matrix saved\n")
cat("  ✓ VSD matrix saved\n", file = summary_file)

# =============================================================================
# Summary
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Summary\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("Original matrix:     %d genes × %d samples\n", 
            nrow(mtx), ncol(mtx)))
cat(sprintf("Filtered matrix:     %d genes × %d samples\n", 
            nrow(filtered_mtx_final_step), ncol(filtered_mtx_final_step)))
cat(sprintf("VSD matrix:          %d genes × %d samples\n", 
            nrow(vsd_matrix), ncol(vsd_matrix)))
cat(sprintf("\nReduction: %.1f%% genes, %.1f%% samples\n",
            (1 - nrow(filtered_mtx_final_step)/nrow(mtx)) * 100,
            (1 - ncol(filtered_mtx_final_step)/ncol(mtx)) * 100))

cat("\nOutput files:\n")
cat(sprintf("  - %s\n", output_filtered))
cat(sprintf("  - %s\n", output_vsd))
cat(sprintf("  - %s\n", output_summary))
cat("\n✓ Processing complete!\n")

# Write summary to file
cat("\n", paste(rep("=", 60), collapse = ""), "\n", file = summary_file)
cat("Summary\n", file = summary_file)
cat(paste(rep("=", 60), collapse = ""), "\n", file = summary_file)
cat(sprintf("Original matrix:     %d genes × %d samples\n", 
            nrow(mtx), ncol(mtx)), file = summary_file)
cat(sprintf("Filtered matrix:     %d genes × %d samples\n", 
            nrow(filtered_mtx_final_step), ncol(filtered_mtx_final_step)), file = summary_file)
cat(sprintf("VSD matrix:          %d genes × %d samples\n", 
            nrow(vsd_matrix), ncol(vsd_matrix)), file = summary_file)
cat(sprintf("\nReduction: %.1f%% genes, %.1f%% samples\n",
            (1 - nrow(filtered_mtx_final_step)/nrow(mtx)) * 100,
            (1 - ncol(filtered_mtx_final_step)/ncol(mtx)) * 100), file = summary_file)

cat("\nOutput files:\n", file = summary_file)
cat(sprintf("  - %s\n", output_filtered), file = summary_file)
cat(sprintf("  - %s\n", output_vsd), file = summary_file)
cat(sprintf("  - %s\n", output_summary), file = summary_file)

close(summary_file)

