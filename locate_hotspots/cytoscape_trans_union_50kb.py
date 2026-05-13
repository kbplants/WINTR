#!/usr/bin/env python3
"""
Build a table of (cluster_id, trans_gene) for Cytoscape from 50 kb winter hotspot regions.

For each region in results/hotspot_groups_winter_50kb.csv, take the union of all
trans genes regulated by the SNPs listed in its `snps` column, using the filtered
trans eQTL file results/winter_trans_eqtl_filtered_p1e-6_min50genes.txt.

Output: results/cytoscape_trans_union_50kb.csv
Columns:
  - cluster_id
  - trans_gene

Each row has exactly one trans_gene; no duplicates per cluster_id.
"""

import os
from collections import defaultdict

import pandas as pd

WORKDIR = "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/024_filter_hotspot"
RESULTS_DIR = os.path.join(WORKDIR, "results")

HOTSPOT_FILE = os.path.join(RESULTS_DIR, "hotspot_groups_winter_50kb.csv")
TRANS_FILE = os.path.join(RESULTS_DIR, "winter_trans_eqtl_filtered_p1e-6_min50genes.txt")
OUT_FILE = os.path.join(RESULTS_DIR, "cytoscape_trans_union_50kb.csv")


def load_snp_genes():
    """Build SNP -> set(trans genes) from filtered trans eQTL file."""
    df = pd.read_csv(TRANS_FILE, sep="\t", usecols=["snps", "gene"])
    snp_genes = defaultdict(set)
    for _, row in df.iterrows():
        snp_genes[row["snps"]].add(row["gene"])
    return snp_genes


def main():
    if not os.path.isfile(HOTSPOT_FILE):
        raise FileNotFoundError(f"Hotspot file not found: {HOTSPOT_FILE}")
    if not os.path.isfile(TRANS_FILE):
        raise FileNotFoundError(f"Trans-eQTL file not found: {TRANS_FILE}")

    os.makedirs(RESULTS_DIR, exist_ok=True)

    snp_genes = load_snp_genes()
    hotspot_df = pd.read_csv(HOTSPOT_FILE)

    rows = []
    for _, row in hotspot_df.iterrows():
        cluster_id = str(row["cluster_id"])
        snp_list = [s.strip() for s in str(row["snps"]).split("|") if s.strip()]
        genes_union = set()
        for snp in snp_list:
            genes_union |= snp_genes.get(snp, set())

        for g in sorted(genes_union):
            rows.append({"cluster_id": cluster_id, "trans_gene": g})

    out_df = pd.DataFrame(rows, columns=["cluster_id", "trans_gene"])
    out_df.to_csv(OUT_FILE, index=False)


if __name__ == "__main__":
    main()

