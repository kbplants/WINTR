# Winter trans eQTL hotspot filter

Filter winter trans eQTL results to high-significance associations and hotspot SNPs (SNPs with many genes). Then identify hotspot groups (nearby SNPs with significant gene overlap) and build region–region overlap matrices, as in `021_locate_hotspot`.

## Working directory

```
/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/024_filter_hotspot
```

## Input

- **Data:** `data/trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new.txt`  
  Winter trans eQTL results (columns: `snps`, `gene`, `statistic`, `pvalue`, `FDR`, `beta`).
- **Fisher test universe:** Total gene count (N) is taken from `021_locate_hotspot/data/winter_expression_matrix.csv` for overlap significance.

## Filters (step 1)

1. **pvalue < 1e-6**  
   Keep only SNP–gene associations with pvalue < 1e-6.

2. **SNPs with ≥ 50 genes**  
   After the pvalue filter, keep only SNPs that have at least 50 unique genes; all associations for those SNPs are retained.

## Scripts (in `scripts/`)

| Script | Description |
|--------|-------------|
| `filter_winter_trans_eqtl.py` | Filter trans eQTL by pvalue and min genes; writes filtered table and summary. |
| `identify_snp_hotspot_groups.py` | Cluster SNPs by distance (50/100/200 kb), keep clusters with significant gene overlap (Fisher p < 0.05); then add one row per SNP that is not in any such hotspot (single-SNP rows) so no SNPs are lost. Writes `hotspot_groups_winter
_*kb.csv` and auxiliary files. |
| `build_hotspot_overlap_matrix.py` | For each distance, build region × region matrix of shared trans gene counts; writes `hotspot_overlap_matrix_winter_*kb.csv`. |
| `cytoscape_trans_union_50kb.py` | From `hotspot_groups_winter_50kb.csv` and the filtered trans eQTL file, build a long table of `(cluster_id, trans_gene)` where trans_gene is the union of all trans genes for SNPs in that region (one gene per row, no duplicates per 
cluster). |
| `cytoscape_type_cis_20kb.py` | For each cluster, find cis genes (from `cis_eqtl_winter_bark_xylem_gene_filter_pf_all_new_e-6.txt`) and 20kb genes (from GFF: genes within 10 kb upstream + 10 kb downstream of each SNP). Output: `cluster_id`, `type` (cis or 20kb), `ge
ne`. |
| `go_enrichment_trans_50kb.R` | GO enrichment on trans genes per cluster from `cytoscape_trans_union_50kb.csv` using clusterProfiler::enricher() with custom annotation from 021_locate_hotspot. Output: `cluster_id`, `go_enrichment` (top 10 GO terms). Run with: `micro
mamba run -n r-env Rscript scripts/go_enrichment_trans_50kb.R` |

**Run order (from working directory):**

```bash
python3 scripts/filter_winter_trans_eqtl.py
python3 scripts/identify_snp_hotspot_groups.py
python3 scripts/build_hotspot_overlap_matrix.py
python3 scripts/cytoscape_trans_union_50kb.py
python3 scripts/cytoscape_type_cis_20kb.py
micromamba run -n r-env Rscript scripts/go_enrichment_trans_50kb.R
```

## Output (results folder)

### From filter

| File | Description |
|------|-------------|
| `winter_trans_eqtl_filtered_p1e-6_min50genes.txt` | Filtered trans eQTL table (same columns as input). |
| `filter_summary.txt` | Counts: total associations, unique SNPs at each step. |

### From hotspot pipeline (same structure as 021_locate_hotspot)

| File | Description |
|------|-------------|
| `hotspot_groups_winter_50kb.csv` | Hotspot regions (50 kb): multi-SNP clusters with significant overlap plus one row per SNP not in any such cluster (n_snps=1, empty overlap stats). Columns: cluster_id, chr, start_bp, end_bp, n_snps, mean_jaccard, mean_shared_trans
_genes, overlap_pvalue, snps. |
| `hotspot_groups_winter_100kb.csv` | Same for 100 kb. |
| `hotspot_groups_winter_200kb.csv` | Same for 200 kb. |
| `hotspot_overlap_matrix_winter_50kb.csv` | Region × region matrix of shared trans gene counts (50 kb). |
| `hotspot_overlap_matrix_winter_100kb.csv` | Same for 100 kb. |
| `hotspot_overlap_matrix_winter_200kb.csv` | Same for 200 kb. |
| `snp_positions_winter_filtered.txt` | SNP ID, chr, pos for filtered SNPs. |
| `snp_trans_gene_counts_winter_filtered.csv` | SNP ID and n_trans_genes. |
| `genomic_clusters_winter_filtered_*kb.txt` | All distance-based clusters (before overlap filter). |
| `cytoscape_trans_union_50kb.csv` | Long-format table for Cytoscape: one row per `(cluster_id, trans_gene)` pair, where trans_gene is from the union of trans genes of all SNPs in that 50 kb region. |
| `cytoscape_type_cis_20kb.csv` | Table with `cluster_id`, `type` (cis or 20kb), `gene`. cis = genes from cis eQTL (pvalue < 1e-6) for SNPs in the cluster; 20kb = genes from GFF within 20 kb of any SNP (10 kb up/downstream). |
| `go_enrichment_trans_50kb.csv` | Table with `cluster_id`, `go_enrichment`. Top 10 GO terms (by adjusted p-value) per cluster from clusterProfiler enrichment on trans genes. |

## Summary (latest run)
- **Unique SNPs remaining (after filter):** **595**
- Associations after pvalue filter: 714,390  
- Associations in final filtered file: 76,916 (only SNPs with ≥ 50 genes)
- **Hotspot groups (overlap p < 0.05):** 77 (50 kb), 79 (100 kb), 73 (200 kb)
- **Single-SNP rows added (so no SNPs lost):** 166 (50 kb), 150 (100 kb), 143 (200 kb)
- **Total rows in table (regions + singles):** 243 (50 kb), 229 (100 kb), 216 (200 kb)
