#!/usr/bin/env python3
"""
Script to rename columns in gene_count_matrix_summer_vsd.csv
- Replace dots (.) with dashes (-) in sample names
- Strip quotes from column names
- Save as new file
"""

import pandas as pd
import sys
from pathlib import Path

# Input and output paths
base_dir = Path("/mnt/winterprojectceph/winter_pipe/combine_batch_results/gf_cnt_mtx")
input_file = base_dir / "results" / "gene_count_matrix_summer_vsd.csv"
output_file = base_dir / "results" / "gene_count_matrix_summer_vsd_renamed.csv"

print("="*60)
print("Renaming columns in summer VSD matrix")
print("="*60)
print(f"\nInput file: {input_file}")
print(f"Output file: {output_file}")

# Read the CSV file
print("\nReading input file...")
df = pd.read_csv(input_file, index_col=0)

print(f"Original matrix: {df.shape[0]} genes × {df.shape[1]} samples")
print(f"\nOriginal column names (first 5):")
print(list(df.columns[:5]))

# Rename columns: replace . with - and strip quotes
print("\nRenaming columns...")
new_columns = []
for col in df.columns:
    # Strip quotes if present
    col_clean = col.strip('"').strip("'")
    # Replace dots with dashes
    col_renamed = col_clean.replace('.', '-')
    new_columns.append(col_renamed)

df.columns = new_columns

print(f"Renamed column names (first 5):")
print(list(df.columns[:5]))

# Save to new file
print(f"\nSaving to: {output_file}")
df.to_csv(output_file)

print(f"\n✓ File saved successfully!")
print(f"  Output: {df.shape[0]} genes × {df.shape[1]} samples")
print(f"  File size: {output_file.stat().st_size / (1024*1024):.1f} MB")




