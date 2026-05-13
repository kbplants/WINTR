# Gene Filtering Pipeline for Summer Count Matrix

## Overview

This directory contains scripts and results for filtering the summer gene count matrix through a comprehensive multi-step filtering pipeline. The pipeline removes low-quality genes and samples, and performs DESeq2 normalization and variance stabilization.

## Directory Structure

```
gf_cnt_mtx/
├── scripts/
│   └── gene_filtering_summer.R    # R script for gene filtering
├── results/
│   ├── gene_count_matrix_summer_filtered.csv  # Filtered count matrix
│   ├── gene_count_matrix_summer_vsd.csv        # VSD normalized matrix
│   └── gene_filtering_summary.txt              # Filtering summary statistics
└── README.md                     # This file
```

## Input File

- **Input**: `/mnt/winterprojectceph/winter_pipe/combine_batch_results/raw_cnt_mtx/results/gene_count_matrix_summer.csv`
  - Format: CSV with gene IDs as first column, samples as remaining columns
  - Contains raw gene counts from RNA-seq analysis

## Gene Filtering Steps

The pipeline performs the following filtering steps in sequence:

### Step 1: Initial Zero Removal
- **Purpose**: Remove genes and samples with no expression
- **Method**: 
  - Replace NA values with 0
  - Remove genes where row sum = 0
  - Remove samples where column sum = 0
- **Rationale**: Eliminates completely unexpressed genes and failed samples

### Step 2: IQR-based Sample Filtering
- **Purpose**: Remove samples with too few expressed genes (outliers)
- **Method**:
  - Calculate number of genes expressed per sample (`colSums(mtx > 0)`)
  - Calculate IQR (Interquartile Range) statistics:
    - Q1 (25th percentile)
    - Q3 (75th percentile)
    - IQR = Q3 - Q1
  - Set threshold: `Q1 - 1.5 × IQR`
  - Keep samples with expressed genes > threshold
- **Rationale**: Removes samples with unusually low gene detection (potential quality issues)

### Step 3: IQR-based Gene Filtering
- **Purpose**: Remove genes expressed in too few samples
- **Method**:
  - Calculate number of samples where each gene is expressed (`rowSums(mtx > 0)`)
  - Calculate IQR statistics for gene expression frequency
  - Set threshold: `Q1 - 1.5 × IQR`
  - Keep genes expressed in more samples than threshold
- **Rationale**: Removes genes that are rarely detected (likely noise or low-quality)

### Step 4: Zero Percentage Filtering (0-10%)
- **Purpose**: Keep only genes with low zero expression across samples
- **Method**:
  - Calculate percentage of zeros per gene: `rowMeans(mtx == 0) × 100`
  - Bin genes into 10% intervals (0-10%, 10-20%, ..., 90-100%)
  - Keep only genes in the 0-10% zero bin
- **Rationale**: Ensures genes are consistently expressed across most samples

### Step 5: 10th Percentile Average Expression Filter
- **Purpose**: Remove genes with very low average expression
- **Method**:
  - Calculate average expression per gene: `rowMeans(mtx)`
  - Calculate 10th percentile threshold
  - Keep genes with average expression > 10th percentile
- **Rationale**: Removes genes with consistently low expression levels

### Step 6: Minimum Expression Filter (>15)
- **Purpose**: Additional filtering based on minimum expression threshold
- **Method**:
  - Load gene statistics from `gene_stats_wint.xlsx` (if available)
  - Filter genes with `min_expr > 15`
  - Keep only genes present in the filtered list
- **Rationale**: Ensures genes have meaningful minimum expression levels
- **Note**: This step is optional and skipped if the gene stats file is not found

### Step 7: DESeq2 Normalization and Variance Stabilization
- **Purpose**: Normalize counts and stabilize variance for downstream analysis
- **Method**:
  - Create DESeqDataSet with design `~ 1` (no covariates)
  - Estimate size factors using `poscounts` method
  - Apply Variance Stabilizing Transformation (VST)
  - Extract VSD matrix
- **Rationale**: 
  - Normalizes for sequencing depth differences
  - Stabilizes variance across expression range
  - Makes data suitable for PCA, clustering, and other downstream analyses

## Output Files

### 1. `gene_count_matrix_summer_filtered.csv`
- **Description**: Filtered gene count matrix after all filtering steps
- **Format**: CSV with `gene_id` as first column
- **Content**: Raw counts for genes and samples that passed all filters
- **Use**: For count-based analyses (e.g., DESeq2 differential expression)

### 2. `gene_count_matrix_summer_vsd.csv`
- **Description**: Variance-stabilized transformed matrix
- **Format**: CSV with `gene_id` as first column
- **Content**: VSD-normalized expression values
- **Use**: For PCA, clustering, correlation analysis, WGCNA

### 3. `gene_filtering_summary.txt`
- **Description**: Detailed summary of filtering statistics
- **Content**: 
  - Number of genes/samples at each step
  - Thresholds used
  - Number of genes/samples removed at each step
  - Final matrix dimensions

## Results

### Generated Files

The pipeline has been executed and generated the following output files:

1. **`gene_count_matrix_summer_filtered.csv`** (56 MB)
   - Filtered count matrix with **21,111 genes × 632 samples**
   - Contains raw counts after all filtering steps
   - File size: 56 MB

2. **`gene_count_matrix_summer_vsd.csv`** (216 MB)
   - VSD normalized matrix with **21,111 genes × 632 samples**
   - Contains variance-stabilized transformed values
   - File size: 216 MB

3. **`gene_filtering_summary.txt`** (2.5 KB)
   - Summary of filtering statistics
   - Contains detailed information about each filtering step

### Filtering Summary

**Original Matrix:**
- **41,133 genes × 641 samples**

**Filtered Matrix:**
- **21,111 genes × 632 samples**

**Reduction:**
- **48.7% genes removed** (20,022 genes)
- **1.4% samples removed** (9 samples)

### Filtering Steps Applied

1. **Initial zero removal**
   - Removed genes and samples with all zeros

2. **IQR-based sample filtering**
   - **Removed 9 samples**
   - Threshold: 25,807 genes expressed
   - Q1: 27,739, Q3: 29,027, IQR: 1,288
   - Samples kept: 632

3. **IQR-based gene filtering**
   - **No genes removed** (threshold was negative: -549)
   - Q1: 165, Q3: 641, IQR: 476
   - All 41,133 genes passed this filter

4. **Zero percentage filter (0-10%)**
   - **Kept 23,457 genes**
   - Removed 17,676 genes with >10% zeros

5. **10th percentile filter**
   - **Kept 21,111 genes**
   - Threshold: 75.41 average expression
   - Removed 2,346 genes below threshold

6. **Min expression filter (>15)**
   - **Skipped** (gene stats file not found: `~/Downloads/gene_stats_wint.xlsx`)
   - Used previous filtered matrix (21,111 genes)

7. **DESeq2 VSD transformation**
   - Applied to filtered matrix
   - Final VSD matrix: **21,111 genes × 632 samples**

### File Locations

All results are saved in:
```
/mnt/winterprojectceph/winter_pipe/combine_batch_results/gf_cnt_mtx/results/
```

## Usage

### From R console:
```r
source("/mnt/winterprojectceph/winter_pipe/combine_batch_results/gf_cnt_mtx/scripts/gene_filtering_summer.R")
```

### From command line:
```bash
Rscript /mnt/winterprojectceph/winter_pipe/combine_batch_results/gf_cnt_mtx/scripts/gene_filtering_summer.R
```

## Requirements

The script requires the following R packages:
- `DESeq2`: Differential expression analysis and normalization
- `dplyr`: Data manipulation
- `readr`: Reading CSV files
- `tidyr`: Data tidying
- `readxl`: Reading Excel files (for optional gene stats)

Install missing packages with:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages(c("dplyr", "readr", "tidyr", "readxl"))
```

## Parameters

The script uses the following default parameters:

- **IQR multiplier**: 1.5 (standard for outlier detection)
- **Zero percentage bin**: 0-10% (genes with ≤10% zeros)
- **Percentile threshold**: 10th percentile (for average expression)
- **Min expression threshold**: 15 (if gene stats file provided)
- **DESeq2 size factor method**: `poscounts` (handles zeros well)

These can be modified in the script if needed.

## Filtering Statistics

After running the script, check `gene_filtering_summary.txt` for:
- Original matrix dimensions
- Dimensions after each filtering step
- Number of genes/samples removed at each step
- Final filtered matrix dimensions
- Percentage reduction in genes and samples

## Notes

1. **Gene Stats File**: The minimum expression filter (Step 6) is optional. If `~/Downloads/gene_stats_wint.xlsx` is not found, the script will skip this step and continue with the previous filtered matrix.

2. **Memory Usage**: Large matrices may require significant memory. Ensure sufficient RAM is available.

3. **Processing Time**: The pipeline may take several minutes to hours depending on matrix size and system resources.

4. **VSD Transformation**: The VSD matrix is recommended for downstream analyses that assume normal distribution (PCA, clustering, correlation).

5. **Count Matrix**: The filtered count matrix should be used for count-based analyses (DESeq2, edgeR, etc.).

## Troubleshooting

### Error: "cannot open file"
- Ensure the input file exists at the specified path
- Check file permissions

### Error: "package 'DESeq2' is not installed"
- Install BiocManager and DESeq2 as shown in Requirements

### Warning: "Gene stats file not found"
- This is expected if the gene stats file is not available
- The script will continue without the min expression filter

### Memory errors
- Reduce matrix size before filtering
- Increase system RAM
- Process in batches

## Author

Generated for winter project analysis pipeline

## Date

2025-01-02


