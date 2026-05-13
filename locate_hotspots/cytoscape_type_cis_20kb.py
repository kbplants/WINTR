#!/usr/bin/env python3
"""
Build cytoscape_type_cis_20kb.csv: cluster_id, type, gene.

For each cluster in hotspot_groups_winter_50kb.csv:
  - cis: genes from cis_eqtl_winter_bark_xylem_gene_filter_pf_all_new_e-6.txt
    (cis eQTL associations for SNPs in the cluster)
  - 20kb: genes from GFF that fall within 20kb of any SNP (10kb upstream + 10kb downstream)

Output: results/cytoscape_type_cis_20kb.csv
"""

import os
import re
from collections import defaultdict

import pandas as pd

WORKDIR = "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/024_filter_hotspot"
DATA_DIR = os.path.join(WORKDIR, "data")
RESULTS_DIR = os.path.join(WORKDIR, "results")

HOTSPOT_FILE = os.path.join(RESULTS_DIR, "hotspot_groups_winter_50kb.csv")
CIS_FILE = os.path.join(DATA_DIR, "cis_eqtl_winter_bark_xylem_gene_filter_pf_all_new_e-6.txt")
GFF_FILE = os.path.join(DATA_DIR, "Ptrichocarpa_444_v3.1.gene_exons.gff3")
OUT_FILE = os.path.join(RESULTS_DIR, "cytoscape_type_cis_20kb.csv")

WINDOW_KB = 10  # 10 kb upstream and downstream = 20 kb total


def parse_snp_id(snp_id):
    """Chr16_13502949 -> (Chr16, 13502949); scaffold_26_488794 -> (scaffold_26, 488794)."""
    s = str(snp_id).strip()
    parts = s.split("_")
    if len(parts) < 2:
        return None, None
    try:
        pos = int(parts[-1])
        chrom = "_".join(parts[:-1])
        return chrom, pos
    except ValueError:
        return None, None


def load_cis_by_snp():
    """SNP -> set of cis gene IDs."""
    df = pd.read_csv(CIS_FILE, sep="\t", usecols=["snps", "gene"])
    snp_genes = defaultdict(set)
    for _, row in df.iterrows():
        snp_genes[row["snps"]].add(row["gene"])
    return dict(snp_genes)


def load_genes_by_chr():
    """Parse GFF: chr -> list of (start, end, gene_id) for type=gene."""
    chr_genes = defaultdict(list)
    with open(GFF_FILE) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, _, ftype, start, end = parts[0], parts[1], parts[2], int(parts[3]), int(parts[4])
            if ftype != "gene":
                continue
            m = re.search(r"ID=([^;]+)", parts[8])
            if m:
                gene_id = m.group(1)
                chr_genes[seqid].append((start, end, gene_id))
    return chr_genes


def genes_within_20kb(chr_genes, chrom, pos):
    """Genes on chrom whose (start,end) overlaps [pos-10kb, pos+10kb]."""
    window_start = max(1, pos - WINDOW_KB * 1000)
    window_end = pos + WINDOW_KB * 1000
    out = set()
    for start, end, gene_id in chr_genes.get(chrom, []):
        if end >= window_start and start <= window_end:
            out.add(gene_id)
    return out


def main():
    if not os.path.isfile(HOTSPOT_FILE):
        raise FileNotFoundError(f"Hotspot file not found: {HOTSPOT_FILE}")
    if not os.path.isfile(CIS_FILE):
        raise FileNotFoundError(f"Cis eQTL file not found: {CIS_FILE}")
    if not os.path.isfile(GFF_FILE):
        raise FileNotFoundError(f"GFF file not found: {GFF_FILE}")

    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("Loading cis eQTL...")
    cis_by_snp = load_cis_by_snp()
    print("Loading GFF genes...")
    chr_genes = load_genes_by_chr()
    print("Loading hotspot groups...")
    hotspot_df = pd.read_csv(HOTSPOT_FILE)

    seen = set()
    rows = []

    for _, row in hotspot_df.iterrows():
        cluster_id = str(row["cluster_id"])
        snp_list = [s.strip() for s in str(row["snps"]).split("|") if s.strip()]

        cis_genes = set()
        for snp in snp_list:
            cis_genes |= cis_by_snp.get(snp, set())

        genes_20kb = set()
        for snp in snp_list:
            chrom, pos = parse_snp_id(snp)
            if chrom is not None and pos is not None:
                genes_20kb |= genes_within_20kb(chr_genes, chrom, pos)

        for g in cis_genes:
            key = (cluster_id, "cis", g)
            if key not in seen:
                seen.add(key)
                rows.append({"cluster_id": cluster_id, "type": "cis", "gene": g})

        for g in genes_20kb:
            key = (cluster_id, "20kb", g)
            if key not in seen:
                seen.add(key)
                rows.append({"cluster_id": cluster_id, "type": "20kb", "gene": g})

    out_df = pd.DataFrame(rows, columns=["cluster_id", "type", "gene"])
    out_df.to_csv(OUT_FILE, index=False)
    print(f"Wrote {OUT_FILE} ({len(out_df)} rows)")


if __name__ == "__main__":
    main()
