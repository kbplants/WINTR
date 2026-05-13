# DESeq2 Analysis Without Expression Filtering

## Overview

This directory contains DESeq2 differential expression analysis results for Summer vs Winter comparison **without applying expression filtering thresholds**. Unlike the standard analysis which filters genes based on:
- Expression in ≥33% of samples
- Standard deviation > 1
- Maximum expression value ≥ 5

This analysis includes **nearly all genes** in the dataset, providing a comprehensive view of differential expression without pre-filtering.

**Note**: A minimal filter is applied to remove genes with zero counts in ALL samples. This is required for DESeq2's size factor estimation and is different from expression-based filtering. Genes with any expression in at least one sample are retained.

## Analysis Details

### Scripts Used

1. **`deseq2_summer_vs_winter_no_threshold.R`**
   - Main DESeq2 analysis script without filtering
   - Location: `/mnt/winterprojectceph/new_winter/001_deseq/test/deseq2_summer_vs_winter_no_threshold.R`
   - Performs DESeq2 analysis on all genes (no expression filtering)

2. **`extract_winter_biased_genes_no_threshold.R`**
   - Extracts winter-biased genes from DESeq2 results
   - Location: `/mnt/winterprojectceph/new_winter/001_deseq/test/extract_winter_biased_genes_no_threshold.R`
   - Filters for FDR < 0.05 and log2FC > 0

### Key Differences from Filtered Analysis

| Aspect | Filtered Analysis | No Threshold Analysis |
|--------|------------------|----------------------|
| Gene filtering | ≥33% samples, sd>1, max≥5 | Minimal: only remove genes with zero counts in ALL samples |
| Total genes tested | ~30,264 | ~42,951 (minus genes with all zeros) |
| Output naming | Standard names | `_no_threshold` suffix |

## Input Files

- **Summer counts**: `/mnt/winterprojectceph/winter_pipe/combine_batch_results/raw_cnt_mtx/results/gene_count_matrix_summer.csv`
- **Winter counts**: `/mnt/winterprojectceph/new_winter/results/raw_count_matrix/gene_count_matrix_winter.csv`
- **Common samples**: 229 samples present in both datasets

## Output Files

All output files are located in this directory (`results/without_threshold/`) and have the `_no_threshold` suffix:

### Main Results Files

1. **`DESeq2results_winter_VS_summer_filtered_FDR0.001_no_threshold.txt`**
   - All genes with FDR < 0.001
   - Columns: gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

2. **`DESeq2results_winter_VS_summer_filtered_FDR0.001_log2FC1_no_threshold.txt`**
   - Genes with FDR < 0.001 and |log2FC| >= 1

3. **`winter_biased_genes_FDR0.05_no_threshold.txt`**
   - Winter-biased genes (FDR < 0.05, log2FC > 0)
   - Extracted using `extract_winter_biased_genes_no_threshold.R`

4. **`winter_biased_genes_FDR0.05_log2FC1_no_threshold.txt`**
   - Strict winter-biased genes (FDR < 0.05, log2FC >= 1)
   - At least 2-fold higher expression in winter

### Count Matrix

- **`gene_count_matrix_summer_winter_common_samples_no_threshold.txt`**
   - Unfiltered count matrix with all genes
   - Used as input for DESeq2

### Diagnostic Plots

- `dispersion_estimates_no_threshold.png` - Dispersion estimates plot
- `MA_plot_all_genes_no_threshold.png` - MA plot for all genes
- `pvalue_distribution_no_threshold.png` - P-value distribution
- `library_sizes_barplot_no_threshold.png` - Library size barplot
- `DESeq2results_winter_VS_summer_MA_plot_no_threshold.png` - MA plot for contrast

### Summary Files

- **`analysis_summary_no_threshold.txt`** - Summary statistics of the analysis
- **`deseq2_no_threshold.log`** - Log file from DESeq2 analysis

## Running the Analysis

### Step 1: Run DESeq2 Analysis (No Filtering)

```bash
cd /mnt/winterprojectceph/new_winter/001_deseq/test
Rscript deseq2_summer_vs_winter_no_threshold.R
```

This will:
- Load summer and winter count matrices
- Find common samples
- Combine datasets (no filtering applied)
- Run DESeq2 analysis on all genes
- Generate results and diagnostic plots

### Step 2: Extract Winter-Biased Genes

```bash
cd /mnt/winterprojectceph/new_winter/001_deseq/test
Rscript extract_winter_biased_genes_no_threshold.R
```

This will:
- Load DESeq2 results from Step 1
- Extract winter-biased genes (FDR < 0.05, log2FC > 0)
- Extract strict winter-biased genes (FDR < 0.05, log2FC >= 1)
- Save results with `_no_threshold` suffix

## Results Summary

The analysis includes nearly all genes from the combined dataset without expression-based pre-filtering. Only genes with zero counts in ALL samples are removed (required for DESeq2). This allows detection of:
- Lowly expressed genes that may be biologically relevant
- Genes with variable expression patterns
- Comprehensive differential expression across the entire transcriptome
- Genes expressed in very few samples (which would be filtered out by the 33% threshold)

### Comparison with Filtered Analysis

| Metric | Filtered | No Threshold |
|--------|----------|--------------|
| Genes tested | ~30,264 | ~42,951 |
| Genes with padj < 0.05 | ~27,995 | (See analysis_summary) |
| Winter-biased (FDR<0.05, log2FC>0) | 9,439 | (See output files) |

## Notes

- **No expression filtering**: All genes are included regardless of expression level
- **DESeq2 internal filtering**: DESeq2 still applies its own independent filtering based on mean counts
- **File naming**: All output files have `_no_threshold` suffix to distinguish from filtered analysis
- **Output directory**: All results are saved in `results/without_threshold/`

## Date

Analysis created: 2025-01-28
