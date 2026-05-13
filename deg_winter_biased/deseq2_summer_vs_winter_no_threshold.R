rm(list=ls())

#############################################################################
####### differential expression analysis of RNAseq data using DESeq2 ########
######## Adapted from R Ployet (ORNL) for Summer vs Winter comparison ########
######## NO FILTERING VERSION - All genes included #########################
#############################################################################

setwd("/mnt/winterprojectceph/new_winter/001_deseq/test")

library("DESeq2")

# Create output directory
output_dir <- "results/without_threshold"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#===============================================================================
# Step 1: Load count matrices
#===============================================================================
cat("Loading count matrices...\n")

# Input files
SUMMER_FILE <- "/mnt/winterprojectceph/winter_pipe/combine_batch_results/raw_cnt_mtx/results/gene_count_matrix_summer.csv"
WINTER_FILE <- "/mnt/winterprojectceph/new_winter/results/raw_count_matrix/gene_count_matrix_winter.csv"

# Load summer counts
summer_counts <- read.csv(SUMMER_FILE, header = TRUE, row.names = 1, check.names = FALSE)
cat("Summer: Loaded", nrow(summer_counts), "genes and", ncol(summer_counts), "samples\n")

# Load winter counts
winter_counts <- read.csv(WINTER_FILE, header = TRUE, row.names = 1, check.names = FALSE)
cat("Winter: Loaded", nrow(winter_counts), "genes and", ncol(winter_counts), "samples\n")

# Convert summer sample names: dots (.) to hyphens (-) to match winter format
colnames(summer_counts) <- gsub("\\.", "-", colnames(summer_counts))
cat("Converted summer sample names: dots to hyphens\n")

#===============================================================================
# Step 2: Find common samples
#===============================================================================
cat("\nIdentifying common samples...\n")

summer_samples <- colnames(summer_counts)
winter_samples <- colnames(winter_counts)
common_samples <- intersect(summer_samples, winter_samples)

cat("Total summer samples:", length(summer_samples), "\n")
cat("Total winter samples:", length(winter_samples), "\n")
cat("Common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("ERROR: No common samples found between summer and winter datasets!")
}

# Filter to common samples only
summer_counts_filtered <- summer_counts[, common_samples, drop = FALSE]
winter_counts_filtered <- winter_counts[, common_samples, drop = FALSE]

cat("Summer filtered:", nrow(summer_counts_filtered), "genes,", ncol(summer_counts_filtered), "samples\n")
cat("Winter filtered:", nrow(winter_counts_filtered), "genes,", ncol(winter_counts_filtered), "samples\n")

#===============================================================================
# Step 3: Combine datasets
#===============================================================================
cat("\nCombining datasets...\n")

# Get all genes (union - keep all genes)
all_genes <- unique(c(rownames(summer_counts_filtered), rownames(winter_counts_filtered)))
cat("Total unique genes:", length(all_genes), "\n")

# Create combined count matrix
combined_counts <- matrix(0, nrow = length(all_genes), ncol = length(common_samples) * 2)
rownames(combined_counts) <- all_genes
colnames(combined_counts) <- c(paste0("summer_", common_samples), paste0("winter_", common_samples))

# Fill in summer counts
summer_genes <- intersect(all_genes, rownames(summer_counts_filtered))
combined_counts[summer_genes, paste0("summer_", common_samples)] <- as.matrix(summer_counts_filtered[summer_genes, common_samples, drop = FALSE])
cat("Added", length(summer_genes), "summer genes\n")

# Fill in winter counts
winter_genes <- intersect(all_genes, rownames(winter_counts_filtered))
combined_counts[winter_genes, paste0("winter_", common_samples)] <- as.matrix(winter_counts_filtered[winter_genes, common_samples, drop = FALSE])
cat("Added", length(winter_genes), "winter genes\n")

# Convert to integer (counts should be integers)
combined_counts <- round(combined_counts)
combined_counts[combined_counts < 0] <- 0
# Replace any NA values with 0
combined_counts[is.na(combined_counts)] <- 0

cat("Combined matrix:", nrow(combined_counts), "genes,", ncol(combined_counts), "samples\n")

#===============================================================================
# Step 4: MINIMAL FILTERING - Remove only genes with zero counts in ALL samples
#===============================================================================
cat("\nApplying minimal filter: removing genes with zero counts in ALL samples...\n")
cat("This is required for DESeq2 size factor estimation.\n")
cat("Total genes before minimal filter:", nrow(combined_counts), "\n")

# Remove genes that have zero counts in ALL samples (required for DESeq2)
# This is different from expression filtering - we're just removing genes with no data
genes_with_expression <- rowSums(combined_counts) > 0
countData <- combined_counts[genes_with_expression, , drop = FALSE]

cat("Genes removed (all zeros):", sum(!genes_with_expression), "\n")
cat("Final count matrix (minimal filter only):", dim(countData)[1], "genes,", dim(countData)[2], "samples\n")
cat("Note: No expression-based filtering applied (33% samples, sd>1, max>=5 thresholds removed)\n")

# Save unfiltered count matrix
write.table(countData, file=file.path(output_dir, "gene_count_matrix_summer_winter_common_samples_no_threshold.txt"), 
            sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
cat("Saved unfiltered count matrix\n")

#===============================================================================
# Step 5: Prepare data for DESeq2
#===============================================================================
cat("\nPreparing data for DESeq2...\n")

# Create experiment labels (summer vs winter)
colData <- DataFrame(treatment=factor(rep(c('summer','winter'), each=length(common_samples))))
rownames(colData) <- colnames(countData)
print(data.frame(colData, sample_id=colnames(countData)))

# Reads number per library
png(file.path(output_dir, "library_sizes_barplot_no_threshold.png"), width = 1200, height = 600)
barplot(colSums(countData)*1e-6, names=colnames(countData), ylab="Library size (millions)", 
        cex.names=0.8, las=2, horiz = FALSE)
dev.off()
cat("Saved library size barplot\n")

#===============================================================================
# Step 6: Create DESeqDataSet and run DESeq2
#===============================================================================
cat("\nCreating DESeqDataSet...\n")

# Creation of dds object
ddsFullCountTable <- DESeqDataSetFromMatrix(countData, colData, design=formula(~treatment))
ddsFullCountTable

# Make sure that 'summer' is the reference level
ddsFullCountTable$treatment <- relevel(ddsFullCountTable$treatment, "summer")

cat("Running DESeq2...\n")
cat("Using library size normalization for size factors (handles zeros better)...\n")

# Calculate size factors using library sizes (total counts per sample)
# This method works even when all genes have zeros
lib_sizes <- colSums(counts(ddsFullCountTable))
cat("Library size range:", range(lib_sizes), "\n")

# Handle zero library sizes by setting them to a small positive value
if (any(lib_sizes == 0)) {
  cat("Warning: Found", sum(lib_sizes == 0), "samples with zero library size. Setting to minimum non-zero value.\n")
  min_nonzero <- min(lib_sizes[lib_sizes > 0])
  lib_sizes[lib_sizes == 0] <- min_nonzero
}

# Calculate size factors
size_factors <- lib_sizes / mean(lib_sizes)
sizeFactors(ddsFullCountTable) <- size_factors

cat("Size factor range:", range(sizeFactors(ddsFullCountTable)), "\n")

# Run DESeq2
dds <- DESeq(ddsFullCountTable)

# DESeq results
res <- results(dds)
summary(res)

#===============================================================================
# Step 7: Diagnostic plots
#===============================================================================
cat("\nCreating diagnostic plots...\n")

alpha <- 0.05

# Plot dispersion estimates
png(file.path(output_dir, "dispersion_estimates_no_threshold.png"), width = 800, height = 600)
plotDispEsts(dds, ymin = 0.01)
dev.off()
cat("Saved dispersion estimates plot\n")

# MA plot
png(file.path(output_dir, "MA_plot_all_genes_no_threshold.png"), width = 800, height = 800)
plotMA(dds, alpha=alpha, main="MA Plot: Winter vs Summer (No Filtering)")
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()
cat("Saved MA plot\n")

# P-value distribution
png(file.path(output_dir, "pvalue_distribution_no_threshold.png"), width = 1200, height = 600)
par(mfrow=c(1,2))
hist(res$pvalue, breaks=50, main="P-value distribution", xlab="p-value")
hist(res$padj, breaks=50, main="Adjusted p-value distribution", xlab="padj")
dev.off()
cat("Saved p-value distribution plots\n")

# Number of differentially expressed genes with FDR < 5%, |log2(foldchange)| >= 1
significants <- res[!is.na(res$padj) & res$padj < alpha & abs(res$log2FoldChange) >= 1 & res$pvalue < alpha,]
cat("Number of significant genes (FDR < 0.05, |log2FC| >= 1, pvalue < 0.05):", dim(significants)[1], "\n")

#===============================================================================
# Step 8: Extract detailed results
#===============================================================================
cat("\nExtracting detailed results...\n")

# Meaning of the columns
mcols(res, use.names=TRUE)

# Identify factor names to be used in results
resultsNames(dds)

# Threshold adjusted pvalue
alpha <- 0.001

# Contrast: winter vs summer
cond1 <- "winter"
cond2 <- "summer"

res2 <- results(dds, contrast=c("treatment", cond1, cond2), pAdjustMethod='fdr')

# Save ALL results first (no FDR filtering) - needed for extracting DEGs with different FDR thresholds
all_results <- as.data.frame(res2)
all_results$gene_id <- rownames(all_results)
write.table(all_results, file=file.path(output_dir, paste0('DESeq2results_', cond1, '_VS_', cond2, '_all_genes_no_threshold.txt')), 
            sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Saved all DESeq2 results (no FDR filter):", nrow(all_results), "genes\n")

# Differentially expressed genes with FDR < alpha
significants2 <- all_results[!is.na(all_results$padj) & (all_results$padj < alpha),]
cat("Significants with FDR <", alpha, ":", dim(significants2)[1], "\n")

# Save all significant results (FDR < alpha)
write.table(significants2, file=file.path(output_dir, paste0('DESeq2results_', cond1, '_VS_', cond2, '_filtered_FDR', alpha, '_no_threshold.txt')), 
            sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Saved results file (FDR <", alpha, ")\n")

# Filter: |log2(foldchange)| >= 1
significants2_fc <- significants2[!is.na(significants2$padj) & (significants2$padj < alpha) & abs(significants2$log2FoldChange) >= 1,]
cat("Significants with FDR <", alpha, "and |log2FC| >= 1:", dim(significants2_fc)[1], "\n")

write.table(significants2_fc, file=file.path(output_dir, paste0('DESeq2results_', cond1, '_VS_', cond2, '_filtered_FDR', alpha, '_log2FC1_no_threshold.txt')), 
            sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Saved results file (FDR <", alpha, ", |log2FC| >= 1)\n")

# Visualize the differentially expressed genes
png(file.path(output_dir, paste0('DESeq2results_', cond1, '_VS_', cond2, '_MA_plot_no_threshold.png')), width = 800, height = 800)
plotMA(res2, alpha=alpha, main=paste(cond1, 'vs', cond2, ' (No Filtering)', sep=' '))
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()
cat("Saved MA plot for contrast\n")

#===============================================================================
# Step 9: Create summary statistics
#===============================================================================
cat("\nCreating summary...\n")

summary_text <- paste0(
  "===============================================================================\n",
  "DESeq2 Analysis: Summer vs Winter (Common Samples) - NO FILTERING VERSION\n",
  "===============================================================================\n",
  "Analysis Date: ", Sys.Date(), "\n\n",
  "Input Files:\n",
  "  Summer counts: ", SUMMER_FILE, "\n",
  "    Genes: ", nrow(summer_counts), "\n",
  "    Samples: ", ncol(summer_counts), "\n\n",
  "  Winter counts: ", WINTER_FILE, "\n",
  "    Genes: ", nrow(winter_counts), "\n",
  "    Samples: ", ncol(winter_counts), "\n\n",
  "Common Samples: ", length(common_samples), "\n\n",
  "Filtering:\n",
  "  MINIMAL FILTER ONLY: Removed genes with zero counts in ALL samples (required for DESeq2)\n",
  "  Total genes before minimal filter: ", dim(combined_counts)[1], "\n",
  "  Genes removed (all zeros): ", sum(!genes_with_expression), "\n",
  "  Genes tested: ", nrow(res2), "\n",
  "  NOTE: No expression-based filtering applied (33% samples, sd>1, max>=5 thresholds removed)\n\n",
  "DESeq2 Results (Winter vs Summer):\n",
  "  Total genes tested: ", nrow(res2), "\n",
  "  Genes with padj < 0.001: ", sum(!is.na(res2$padj) & res2$padj < 0.001), "\n",
  "  Genes with padj < 0.01: ", sum(!is.na(res2$padj) & res2$padj < 0.01), "\n",
  "  Genes with padj < 0.05: ", sum(!is.na(res2$padj) & res2$padj < 0.05), "\n",
  "  Genes with padj < 0.001 and |log2FC| >= 1: ", dim(significants2_fc)[1], "\n\n",
  "Output Files:\n",
  "  Unfiltered count matrix: ", file.path(output_dir, "gene_count_matrix_summer_winter_common_samples_no_threshold.txt"), "\n",
  "  Results (FDR < 0.001): ", file.path(output_dir, paste0("DESeq2results_winter_VS_summer_filtered_FDR0.001_no_threshold.txt")), "\n",
  "  Results (FDR < 0.001, |log2FC| >= 1): ", file.path(output_dir, paste0("DESeq2results_winter_VS_summer_filtered_FDR0.001_log2FC1_no_threshold.txt")), "\n",
  "  Diagnostic plots: ", file.path(output_dir, "dispersion_estimates_no_threshold.png"), ", ", file.path(output_dir, "MA_plot_all_genes_no_threshold.png"), ", ", file.path(output_dir, "pvalue_distribution_no_threshold.png"), "\n",
  "  Library sizes: ", file.path(output_dir, "library_sizes_barplot_no_threshold.png"), "\n",
  "  MA plot for contrast: ", file.path(output_dir, paste0("DESeq2results_winter_VS_summer_MA_plot_no_threshold.png")), "\n\n",
  "===============================================================================\n"
)

writeLines(summary_text, file.path(output_dir, "analysis_summary_no_threshold.txt"))
cat(summary_text)

cat("\n===============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("===============================================================================\n")
cat("Analysis completed:", Sys.time(), "\n")
