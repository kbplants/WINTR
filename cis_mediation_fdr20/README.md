# 032 cis mediation FDR 0.2 filtering

This directory contains mediation / eQTL results and helper scripts.

## Data

- **Input files**: Located in the `data` directory, e.g.  
  - `cis_eqtl_summer_xylem_gene_filter_pf_all_new.txt`  
  - `cis_eqtl_winter_bark_xylem_gene_filter_pf_all_new.txt`  
  - `trans_eqtl_results_summer_xylem_gene_filter_pf_all_new.txt`  
  - `trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new.txt`

  These are tab-delimited text files with a header row that includes an `FDR` column.

- **FDR 0.2-filtered outputs**: Also in `data`, with the same base name plus an `_fdr_20` suffix, e.g.  
  - `cis_eqtl_summer_xylem_gene_filter_pf_all_new_fdr_20.txt`

  Each `_fdr_20` file contains only rows where `FDR <= 0.2`, and preserves the original header.

- **SNPs with at least 50 trans genes** (trans-eQTL only):  
  For the trans-eQTL result files, there are additional outputs with the suffix `_atleast50genes`, e.g.  
  - `trans_eqtl_results_summer_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt`  
  - `trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt`

  These `_atleast50genes` files are derived from the corresponding `_fdr_20` trans-eQTL files and contain only SNPs that have at least 50 associated (trans) genes (counted as unique `gene` values per SNP).

## Cis–trans mediation (winter)

Mediation analysis for **winter** uses the winter eQTL and expression inputs. Outputs are written under `results/winter/`.

**Inputs (in `data/`):**

- `cis_eqtl_winter_bark_xylem_gene_filter_pf_all_new_fdr_20.txt` — cis-eQTL (FDR ≤ 0.2)
- `trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt` — trans-eQTL (FDR ≤ 0.2, SNPs with ≥50 trans genes)
- `winter_expression_matrix.csv` — expression matrix (must have `gene_id` column)
- `genotype_matrix_output1.txt` — genotype matrix (first column = SNP id, e.g. `snps`)

**Outputs (in `results/winter/`):**

- `cis_trans_triplets_all.tsv` — all SNP–cis gene–trans gene triplets (from `build_triplets.R`)
- `mediation_results.tsv` — mediation regression results (from `run_mediation.R`)

### Running the winter mediation pipeline

From this directory (`032_cis_mediation_fdr20`):

```bash
# 1. Build winter triplets (cis × trans pairs per overlapping SNP)
Rscript scripts/build_triplets.R .

# 2. Run winter mediation regressions and write results
Rscript scripts/run_mediation.R .
```

You can pass an absolute path instead of `.` if you run from elsewhere.

### Winter 20 kb window mediation

This variant uses **genes within ±10 kb (20 kb window)** of each SNP from the winter trans-eQTL file (≥50 trans genes, FDR ≤ 0.2), based on genomic coordinates in the GFF3 file.

**Additional inputs (in `data/`):**

- `trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt` — winter trans-eQTL (as above)
- `Ptrichocarpa_444_v3.1.gene_exons.gff3` — gene annotations (used to find genes within ±10 kb of each SNP)
- **Winter 20 kb mediation** uses the expression matrix in `data/wint_raw_cnt_mtx/gene_count_matrix_winter.csv` (not `winter_expression_matrix.csv`).

**Additional outputs (in `results/winter/`):**

- `cis_trans_triplets_all_20kb_winter.tsv` — SNP–window-gene–trans triplets where the **window gene** lies within 20 kb (±10 kb) of the SNP.
- `mediation_results_20kb.tsv` — mediation regression results for these 20 kb window triplets.

**Running the winter 20 kb pipeline:**
```bash
# 1. Build winter 20 kb SNP–cis-window–trans triplets
Rscript scripts/build_triplets_20kb_winter.R .

# 2. Run winter 20 kb mediation regressions and write results
Rscript scripts/run_mediation_20kb_winter.R .
```

## Cis–trans mediation (summer)

Mediation analysis for **summer** uses the summer eQTL and expression inputs. Outputs are written under `results/summer/`.

**Inputs (in `data/`):**

- `cis_eqtl_summer_xylem_gene_filter_pf_all_new_fdr_20.txt` — cis-eQTL (FDR ≤ 0.2)
- `trans_eqtl_results_summer_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt` — trans-eQTL (FDR ≤ 0.2, SNPs with ≥50 trans genes)
- `summer_expression_data.csv` — expression matrix (must have `gene_id` column)
- `genotype_matrix_output1.txt` — genotype matrix (first column = SNP id, e.g. `snps`)

**Outputs (in `results/summer/`):**

- `cis_trans_triplets_all.tsv` — all SNP–cis gene–trans gene triplets (from `build_triplets_summer.R`)
- `mediation_results.tsv` — mediation regression results (from `run_mediation_summer.R`)

### Running the summer mediation pipeline

From this directory (`032_cis_mediation_fdr20`):

```bash
# 1. Build summer triplets (cis × trans pairs per overlapping SNP)
Rscript scripts/build_triplets_summer.R .

# 2. Run summer mediation regressions and write results
Rscript scripts/run_mediation_summer.R .
```

### Summer 20 kb window mediation

This variant uses **genes within ±10 kb (20 kb window)** of each SNP from the summer trans-eQTL file (≥50 trans genes, FDR ≤ 0.2), based on genomic coordinates in the GFF3 file. The mediator gene is labelled **window_gene** (not cis_gene).

**Additional inputs (in `data/`):**

- `trans_eqtl_results_summer_xylem_gene_filter_pf_all_new_fdr_20_atleast50genes.txt` — summer trans-eQTL (as above)
- `Ptrichocarpa_444_v3.1.gene_exons.gff3` — gene annotations (used to find genes within ±10 kb of each SNP)
- **Summer 20 kb mediation** uses the expression matrix in `data/summ_raw_cnt_mtx/gene_count_matrix_summer.csv` (not `summer_expression_data.csv`).

**Additional outputs (in `results/summer/`):**

- `cis_trans_triplets_all_20kb_summer.tsv` — SNP–window-gene–trans triplets where the **window gene** lies within 20 kb (±10 kb) of the SNP.
- `mediation_results_20kb.tsv` — mediation regression results for these 20 kb window triplets (columns include `window_gene`).

**Running the summer 20 kb pipeline:**

```bash
# 1. Build summer 20 kb SNP–window-gene–trans triplets
Rscript scripts/build_triplets_20kb_summer.R .

# 2. Run summer 20 kb mediation regressions and write results
Rscript scripts/run_mediation_20kb_summer.R .
```


---

## Scripts

- **`scripts/build_triplets.R`**  
  Builds SNP–cis gene–trans gene triplets from the winter cis and trans eQTL files (overlapping SNPs only). Writes `results/winter/cis_trans_triplets_all.tsv`.

- **`scripts/run_mediation.R`**  
  Runs mediation regressions (Sobel-style) for each winter triplet using the winter expression and genotype matrices. Writes `results/winter/mediation_results.tsv`. Requires `data.table` in R.

- **`scripts/build_triplets_20kb_winter.R`**  
  Builds SNP–window-gene–trans triplets from the winter trans-eQTL file (≥50 genes) and GFF3: for each SNP, finds genes within ±10 kb, then joins with trans genes. Writes `results/winter/cis_trans_triplets_all_20kb_winter.tsv` (columns: `snps`, `window_gene`, `trans_
gene`).

- **`scripts/run_mediation_20kb_winter.R`**  
  Runs mediation regressions on the winter 20 kb triplets using the expression matrix in `data/wint_raw_cnt_mtx/gene_count_matrix_winter.csv` and the genotype matrix. Writes `results/winter/mediation_results_20kb.tsv` (columns include `window_gene`). Requires `data.t
able` in R.

- **`scripts/build_triplets_summer.R`**  
  Builds SNP–cis gene–trans gene triplets from the summer cis and trans eQTL files (overlapping SNPs only). Writes `results/summer/cis_trans_triplets_all.tsv`.

- **`scripts/run_mediation_summer.R`**  
  Runs mediation regressions (Sobel-style) for each summer triplet using the summer expression and genotype matrices. Writes `results/summer/mediation_results.tsv`. Requires `data.table` in R.

- **`scripts/build_triplets_20kb_summer.R`**  
  Builds SNP–window-gene–trans triplets from the summer trans-eQTL file (≥50 genes) and GFF3: for each SNP, finds genes within ±10 kb, then joins with trans genes. Writes `results/summer/cis_trans_triplets_all_20kb_summer.tsv` (columns: `snps`, `window_gene`, `trans_
gene`).

- **`scripts/run_mediation_20kb_summer.R`**  
  Runs mediation regressions on the summer 20 kb triplets using the expression matrix in `data/summ_raw_cnt_mtx/gene_count_matrix_summer.csv` and the genotype matrix. Writes `results/summer/mediation_results_20kb.tsv` (columns include `window_gene`). Requires `data.t
able` in R.

- **`scripts/filter_fdr_0_2.py`**  
  Python script that:
  - Filters mediation / eQTL result files to `FDR <= 0.2`.
  - For trans-eQTL result files (those whose names start with `trans_eqtl_results_`), it also generates `_atleast50genes` files that retain only SNPs with at least 50 associated trans genes.

  - **Default behavior** (no arguments):  
    - Looks in `../data` for all `.txt` files.  
    - For each file with an `FDR` column, writes a filtered version with `_fdr_20` inserted before the extension in the same `data` directory.  
    - For each trans-eQTL `_fdr_20` file, writes an additional `_atleast50genes` file that keeps only SNPs with at least 50 trans genes.

  - **Optional usage with explicit files**:  
    - You can pass one or more files to filter:
      ```bash
      python scripts/filter_fdr_0_2.py data/cis_eqtl_summer_xylem_gene_filter_pf_all_new.txt
      ```

### Running the FDR 0.2 filter

From this directory:

```bash
python scripts/filter_fdr_0_2.py
```
