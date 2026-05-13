#!/usr/bin/env Rscript
# =============================================================================
# Script: rename_matrix_columns_to_sample_codes.R
# Purpose: Replace matrix column names (file names) with sample_code names
#          Based on metadata matching from winter_samples.csv
# Author: Generated for winter project
# Date: 2025-01-02
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(data.table)
  library(tibble)
})

# =============================================================================
# Configuration
# =============================================================================

# Input files
gene_count_file <- "/mnt/winterprojectceph/new_winter/results/raw_count_matrix/gene_count_matrix_winter.csv"
samples_file <- "/mnt/winterprojectceph/winter_samples.csv"
metadata_file <- "/mnt/winterprojectceph/winter_pipe/combine_batch_results/raw_cnt_mtx/spReport_query_504713_7_25.2025.csv"

# Output file
output_file <- "/mnt/winterprojectceph/new_winter/results/raw_count_matrix/gene_count_matrix_winter.csv"

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Replace Matrix Column Names with Sample Codes\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# =============================================================================
# Step 1: Load and prepare metadata
# =============================================================================

cat("Step 1: Loading metadata files...\n")

# Load sample file names
winter_samples <- read_csv(samples_file, col_names = FALSE, show_col_types = FALSE)
colnames(winter_samples) <- c("file_name")
cat(sprintf("  Loaded %d sample file names\n", nrow(winter_samples)))

# Load metadata
met_data <- read_csv(metadata_file, show_col_types = FALSE)
cat(sprintf("  Loaded metadata for %d entries\n", nrow(met_data)))

# Transform RQC Seq Unit Name to match file names
# Change .fastq.gz in met_data to .filter-RNA.fastq.gz
met_data <- met_data %>%
  mutate(`RQC Seq Unit Name` = str_replace(`RQC Seq Unit Name`, "\\.fastq\\.gz$", ".filter-RNA.fastq.gz"))

# =============================================================================
# Step 2: Create mapping from file names to sample codes
# =============================================================================

cat("\nStep 2: Creating file name to sample code mapping...\n")

# Match elements in file_name column from winter_samples and RQC Seq Unit Name from met_data
winter_samples_updated <- winter_samples %>%
  left_join(
    met_data %>% select(
      `RQC Seq Unit Name`, 
      `Sample Name`, 
      `Sample Plate/Tube Name`, 
      `Sample Isolated From`,
      `Sample QC Status`
    ), 
    by = c("file_name" = "RQC Seq Unit Name")
  )

# Add season-tissue classification
winter_samples_updated <- winter_samples_updated %>%
  mutate(season_tissue = case_when(
    `Sample Plate/Tube Name` == "WINTR_P1" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P_2" & `Sample Isolated From` == "Bark" ~ "summer_bark",
    `Sample Plate/Tube Name` == "WINTR_P_2" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_P3" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P4" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P5" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P6" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P7" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P8" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P9" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P10" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P11" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_P12" & `Sample Isolated From` == "Bark/xylem" ~ "winter_xylem_bark",
    `Sample Plate/Tube Name` == "WINTR_sum_P1" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_sum_P2" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_sum_P3" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_sum_P4" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_sum_P5" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_sum_P6" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_sum_P7" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_sum_P8" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    `Sample Plate/Tube Name` == "WINTR_sum_P9" & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    str_detect(`Sample Plate/Tube Name`, "WINTR_sum") & `Sample Isolated From` == "xylem" ~ "summer_xylem",
    TRUE ~ NA_character_
  )) %>%
  mutate(season = case_when(
    str_starts(season_tissue, "summer") ~ "summer",
    str_starts(season_tissue, "winter") ~ "winter",
    TRUE ~ NA_character_
  ))

# Replace 'Bark' with 'bark' in 'Sample Isolated From' column
winter_samples_updated <- winter_samples_updated %>%
  mutate(`Sample Isolated From` = str_replace_all(`Sample Isolated From`, "Bark", "bark"))

# Create 'file' column to match matrix column names
# Matrix columns have format: XXXX.filter-RNA.fastq (no .gz, no X prefix)
# Simply remove .gz from file_name to match matrix column names
winter_samples_updated <- winter_samples_updated %>%
  mutate(
    file = str_remove(file_name, "\\.gz$")
  )

# Extract sample_code from Sample Name
winter_samples_updated <- winter_samples_updated %>%
  mutate(sample_code = case_when(
    # Case 1: Sample Name contains "summer" or "winter", and Plate/Tube does NOT contain "SUM"
    str_detect(`Sample Name`, "summer|winter") & 
      !str_detect(`Sample Plate/Tube Name`, regex("sum", ignore_case = TRUE)) ~ 
        str_match(`Sample Name`, "^(?:[^_]+_){3}([^_]+)_")[,2],
    # Case 2: Plate/Tube Name contains "sum" (case-insensitive)
    str_detect(`Sample Plate/Tube Name`, regex("sum", ignore_case = TRUE)) ~ 
        str_match(`Sample Name`, "^(?:[^_]+_){3}(.*)")[,2],
    # Default/fallback
    TRUE ~ str_extract(`Sample Name`, "[^_]+$")
  )) %>%
  mutate(sample_code = str_remove(sample_code, "_S$"))

# Filter out rows without season
winter_samples_updated <- winter_samples_updated %>% 
  filter(!is.na(season))

cat(sprintf("  Created mapping for %d samples\n", nrow(winter_samples_updated)))
cat(sprintf("  Winter samples: %d\n", sum(winter_samples_updated$season == "winter")))
cat(sprintf("  Summer samples: %d\n", sum(winter_samples_updated$season == "summer")))

# =============================================================================
# Step 3: Load gene count matrix
# =============================================================================

cat("\nStep 3: Loading gene count matrix...\n")
cat(sprintf("  Reading: %s\n", gene_count_file))

# Read gene count matrix
gene_count_mtx <- fread(gene_count_file, data.table = FALSE)
rownames(gene_count_mtx) <- gene_count_mtx[, 1]  # Set first column as row names
gene_count_mtx <- gene_count_mtx[, -1, drop = FALSE]  # Remove first column

cat(sprintf("  Matrix dimensions: %d genes × %d samples\n", nrow(gene_count_mtx), ncol(gene_count_mtx)))

# =============================================================================
# Step 4: Map column names to sample codes
# =============================================================================

cat("\nStep 4: Mapping column names to sample codes...\n")

# Get current column names from matrix
current_cols <- colnames(gene_count_mtx)

# Filter metadata to only include files that are in the matrix
winter_samples_filtered <- winter_samples_updated %>%
  filter(file %in% current_cols)

cat(sprintf("  Files in metadata matching matrix columns: %d\n", nrow(winter_samples_filtered)))

# Create mapping from file (matrix column name) to sample_code
file_to_sample_code <- setNames(
  winter_samples_filtered$sample_code, 
  winter_samples_filtered$file
)

# Replace column names in matrix
new_colnames <- colnames(gene_count_mtx)
mapped_count <- 0

for (i in seq_along(new_colnames)) {
  if (new_colnames[i] %in% names(file_to_sample_code)) {
    new_colnames[i] <- file_to_sample_code[new_colnames[i]]
    mapped_count <- mapped_count + 1
  }
}

colnames(gene_count_mtx) <- new_colnames

cat(sprintf("  Mapped %d out of %d columns to sample codes\n", mapped_count, length(current_cols)))

if (mapped_count < length(current_cols)) {
  unmapped <- length(current_cols) - mapped_count
  cat(sprintf("  WARNING: %d columns could not be mapped (kept original file names)\n", unmapped))
}

# =============================================================================
# Step 5: Handle duplicate column names
# =============================================================================

cat("\nStep 5: Handling duplicate column names...\n")

# Check for duplicates
col_counts <- table(colnames(gene_count_mtx))
duplicates <- names(col_counts)[col_counts > 1]

if (length(duplicates) > 0) {
  cat(sprintf("  Found %d duplicate column names\n", length(duplicates)))
  
  # Append suffixes to duplicates
  new_colnames <- colnames(gene_count_mtx)
  counter <- list()
  
  for (i in seq_along(new_colnames)) {
    col_name <- new_colnames[i]
    
    if (col_name %in% duplicates) {
      if (is.null(counter[[col_name]])) {
        counter[[col_name]] <- 1
      } else {
        counter[[col_name]] <- counter[[col_name]] + 1
        new_colnames[i] <- paste0(col_name, "_", counter[[col_name]])
      }
    }
  }
  
  colnames(gene_count_mtx) <- new_colnames
  fixed_count <- sum(colnames(gene_count_mtx) != current_cols)
  cat(sprintf("  Fixed duplicate columns by appending suffixes\n"))
} else {
  cat("  No duplicate column names found\n")
}

# =============================================================================
# Step 6: Save results
# =============================================================================

cat("\nStep 6: Saving results...\n")

# Convert back to data frame with gene_id column
gene_count_mtx_out <- as.data.frame(gene_count_mtx)
gene_count_mtx_out <- cbind(data.frame(gene_id = rownames(gene_count_mtx_out)), gene_count_mtx_out)

cat(sprintf("  Writing to: %s\n", output_file))
fwrite(gene_count_mtx_out, output_file, row.names = FALSE)

cat("  ✓ Matrix saved with sample code column names\n")

# =============================================================================
# Summary
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Summary\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("Original columns: %d\n", length(current_cols)))
cat(sprintf("Columns mapped to sample codes: %d\n", mapped_count))
cat(sprintf("Final matrix: %d genes × %d samples\n", nrow(gene_count_mtx), ncol(gene_count_mtx)))
cat(sprintf("\nOutput file: %s\n", output_file))
cat("\n✓ Processing complete!\n")


