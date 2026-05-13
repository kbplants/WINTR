#!/usr/bin/env python3
"""
Build hotspot region × hotspot region matrices of overlapping trans gene counts
for the filtered winter trans eQTL dataset.

Reads: results/hotspot_groups_winter_{50,100,200}kb.csv
       results/winter_trans_eqtl_filtered_p1e-6_min50genes.txt
For each distance window: rows and columns = hotspot cluster_id; value = number of
trans genes shared between the two regions (union of genes in region i ∩ union of genes in region j).
Saves: results/hotspot_overlap_matrix_winter_{50,100,200}kb.csv
"""

import os
import sys
from collections import defaultdict

import pandas as pd
import numpy as np

WORKDIR = "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/024_filter_hotspot"
os.chdir(WORKDIR)

RESULTS_DIR = "results"
TRANS_FILE = os.path.join(WORKDIR, RESULTS_DIR, "winter_trans_eqtl_filtered_p1e-6_min50genes.txt")
DISTANCE_KB = [50, 100, 200]


def load_snp_genes():
    """SNP -> set of trans gene IDs from filtered trans eQTL."""
    print("Loading filtered trans eQTL (pvalue < 1e-6, >= 50 genes)...")
    df = pd.read_csv(TRANS_FILE, sep="\t", usecols=["snps", "gene"])
    snp_genes = defaultdict(set)
    for _, row in df.iterrows():
        snp_genes[row["snps"]].add(row["gene"])
    print(f"  Loaded {len(snp_genes)} SNPs.")
    return dict(snp_genes)


def genes_for_region(snp_genes, snp_list):
    """Union of trans genes over all SNPs in the region."""
    out = set()
    for s in snp_list:
        out |= snp_genes.get(s, set())
    return out


def main():
    if not os.path.isfile(TRANS_FILE):
        print(f"ERROR: File not found: {TRANS_FILE}", file=sys.stderr)
        sys.exit(1)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    snp_genes = load_snp_genes()

    for distance_kb in DISTANCE_KB:
        hotspot_path = os.path.join(RESULTS_DIR, f"hotspot_groups_winter_{distance_kb}kb.csv")
        if not os.path.isfile(hotspot_path):
            print(f"  Skip {distance_kb} kb: {hotspot_path} not found.", file=sys.stderr)
            continue
        print(f"\n--- Distance window: {distance_kb} kb ---")
        df = pd.read_csv(hotspot_path)
        cluster_ids = df["cluster_id"].astype(str).tolist()
        n = len(cluster_ids)
        print(f"  Hotspot regions: {n}")

        region_genes = {}
        for _, row in df.iterrows():
            cid = str(row["cluster_id"])
            snps = [s.strip() for s in str(row["snps"]).split("|")]
            region_genes[cid] = genes_for_region(snp_genes, snps)

        matrix = np.zeros((n, n), dtype=int)
        for i in range(n):
            gi = region_genes[cluster_ids[i]]
            matrix[i, i] = len(gi)
            for j in range(i + 1, n):
                overlap = len(gi & region_genes[cluster_ids[j]])
                matrix[i, j] = overlap
                matrix[j, i] = overlap

        out_df = pd.DataFrame(matrix, index=cluster_ids, columns=cluster_ids)
        out_path = os.path.join(RESULTS_DIR, f"hotspot_overlap_matrix_winter_{distance_kb}kb.csv")
        out_df.to_csv(out_path)
        print(f"  Wrote {out_path} ({n}x{n})")

    print("\nDone. Overlap matrices in", RESULTS_DIR)


if __name__ == "__main__":
    main()
