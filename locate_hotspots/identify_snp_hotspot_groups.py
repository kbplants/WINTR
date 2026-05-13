#!/usr/bin/env python3
"""
Identify groups of nearby SNPs that regulate the same trans genes (hotspot groups)
for the filtered winter trans eQTL dataset (pvalue < 1e-6, >= 50 genes).

Uses: results/winter_trans_eqtl_filtered_p1e-6_min50genes.txt
Hotspot filter: gene overlap between SNP pairs must be significant (Fisher's exact p < 0.05).
Outputs to results/
"""

import os
import sys
from collections import defaultdict

import pandas as pd
from scipy import stats

WORKDIR = "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/024_filter_hotspot"
os.chdir(WORKDIR)

RESULTS_DIR = "results"
TRANS_FILE = os.path.join(WORKDIR, RESULTS_DIR, "winter_trans_eqtl_filtered_p1e-6_min50genes.txt")
# Use expression matrix from 021 for Fisher's test universe (same as original pipeline)
EXPRESSION_FILE = "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/021_locate_hotspot/data/winter_expression_matrix.csv"

DISTANCE_KB = [50, 100, 200]
MIN_CLUSTER_SIZE = 2
MAX_OVERLAP_PVALUE = 0.05


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


def load_total_genes():
    """Total genes (N) from winter expression matrix (universe for Fisher's test)."""
    df = pd.read_csv(EXPRESSION_FILE, usecols=[0])
    return len(df)


def load_snp_positions_and_genes():
    """Load filtered trans eQTL: derive chr/pos from SNP IDs, build SNP->genes."""
    print("Loading total genes (N) from winter_expression_matrix.csv...")
    total_genes = load_total_genes()
    print(f"  N = {total_genes} genes (universe for Fisher's test)")

    print("Loading filtered trans eQTL (pvalue < 1e-6, >= 50 genes)...")
    df = pd.read_csv(TRANS_FILE, sep="\t", usecols=["snps", "gene"])
    snp_genes = defaultdict(set)
    snp_chr_pos = {}

    for _, row in df.iterrows():
        snp = row["snps"]
        g = row["gene"]
        snp_genes[snp].add(g)
        if snp not in snp_chr_pos:
            chr_name, pos = parse_snp_id(snp)
            if chr_name is not None and pos is not None:
                snp_chr_pos[snp] = (chr_name, pos)

    rows = [{"snp_id": s, "chr": snp_chr_pos[s][0], "pos": snp_chr_pos[s][1]}
            for s in snp_chr_pos if snp_chr_pos[s][0] and snp_chr_pos[s][1] is not None]
    snp_df = pd.DataFrame(rows)
    snp_df = snp_df.drop_duplicates(subset=["snp_id"])
    snp_df = snp_df.sort_values(["chr", "pos"]).reset_index(drop=True)
    print(f"  Found {len(snp_df)} unique SNPs.")
    return snp_df, dict(snp_genes), total_genes


def union_find_cluster(positions, max_dist_bp):
    if len(positions) == 0:
        return []
    positions = sorted(positions)
    parent = {p: p for p in positions}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i, p in enumerate(positions):
        for j in range(i + 1, len(positions)):
            q = positions[j]
            if q - p > max_dist_bp:
                break
            union(p, q)

    clusters = defaultdict(list)
    for p in positions:
        clusters[find(p)].append(p)
    return list(clusters.values())


def cluster_snps_by_distance(snp_df, distance_kb):
    max_dist_bp = distance_kb * 1000
    all_clusters = []
    for chr_name, grp in snp_df.groupby("chr", sort=False):
        pos_to_snp = dict(zip(grp["pos"], grp["snp_id"]))
        positions = grp["pos"].tolist()
        pos_clusters = union_find_cluster(positions, max_dist_bp)
        for pos_list in pos_clusters:
            snp_list = [pos_to_snp[p] for p in pos_list]
            if len(snp_list) >= MIN_CLUSTER_SIZE:
                all_clusters.append({"chr": chr_name, "snps": snp_list, "positions": sorted(pos_list)})
    return all_clusters


def jaccard(a, b):
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def mean_pairwise_jaccard(snp_genes, snp_list):
    sets = [snp_genes.get(s, set()) for s in snp_list]
    n = len(snp_list)
    if n < 2:
        return 0.0
    total = 0.0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            total += jaccard(sets[i], sets[j])
            count += 1
    return total / count if count else 0.0


def mean_pairwise_shared(snp_genes, snp_list):
    sets = [snp_genes.get(s, set()) for s in snp_list]
    n = len(snp_list)
    if n < 2:
        return 0
    total = 0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            total += len(sets[i] & sets[j])
            count += 1
    return total / count if count else 0


def fisher_overlap_pvalue(a_set, b_set, total_genes):
    """Fisher's exact test: is the overlap of A and B significant given total genes?"""
    a, b = len(a_set), len(b_set)
    k = len(a_set & b_set)
    if k == 0:
        return 1.0
    table = [[k, a - k], [b - k, total_genes - a - b + k]]
    if table[1][1] < 0:
        return 1.0
    _, p = stats.fisher_exact(table, alternative="greater")
    return p


def min_overlap_pvalue(snp_genes, snp_list, total_genes):
    """Minimum Fisher p-value across all SNP pairs in cluster."""
    sets = [snp_genes.get(s, set()) for s in snp_list]
    n = len(snp_list)
    if n < 2:
        return None
    min_p = 1.0
    for i in range(n):
        for j in range(i + 1, n):
            p = fisher_overlap_pvalue(sets[i], sets[j], total_genes)
            min_p = min(min_p, p)
    return min_p


def main():
    if not os.path.isfile(TRANS_FILE):
        print(f"ERROR: File not found: {TRANS_FILE}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(EXPRESSION_FILE):
        print(f"ERROR: Expression file not found: {EXPRESSION_FILE}", file=sys.stderr)
        sys.exit(1)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    snp_df, snp_genes, total_genes = load_snp_positions_and_genes()
    snp_df.to_csv(os.path.join(RESULTS_DIR, "snp_positions_winter_filtered.txt"), sep="\t", index=False)

    snp_gene_counts = [(s, len(snp_genes.get(s, set()))) for s in snp_df["snp_id"]]
    pd.DataFrame(snp_gene_counts, columns=["snp_id", "n_trans_genes"]).to_csv(
        os.path.join(RESULTS_DIR, "snp_trans_gene_counts_winter_filtered.csv"), index=False
    )

    for distance_kb in DISTANCE_KB:
        print(f"\n--- Distance window: {distance_kb} kb ---")
        clusters = cluster_snps_by_distance(snp_df, distance_kb)
        print(f"  Genomic clusters (size >= {MIN_CLUSTER_SIZE}): {len(clusters)}")

        rows = []
        for c in clusters:
            snps = c["snps"]
            mean_j = mean_pairwise_jaccard(snp_genes, snps)
            mean_shared = mean_pairwise_shared(snp_genes, snps)
            overlap_pv = min_overlap_pvalue(snp_genes, snps, total_genes)
            chr_name = c["chr"]
            start, end = min(c["positions"]), max(c["positions"])
            rows.append({
                "cluster_id": f"{chr_name}_{start}_{end}",
                "chr": chr_name,
                "start_bp": start,
                "end_bp": end,
                "n_snps": len(snps),
                "mean_jaccard": round(mean_j, 4),
                "mean_shared_trans_genes": round(mean_shared, 1),
                "overlap_pvalue": overlap_pv if overlap_pv is not None else "",
                "snps": "|".join(snps),
            })

        clusters_df = pd.DataFrame(rows)
        clusters_df.to_csv(
            os.path.join(RESULTS_DIR, f"genomic_clusters_winter_filtered_{distance_kb}kb.txt"),
            sep="\t",
            index=False,
        )

        hotspot = clusters_df[
            clusters_df["overlap_pvalue"].notna() &
            (clusters_df["overlap_pvalue"] < MAX_OVERLAP_PVALUE)
        ].copy()

        # SNPs that appear in at least one multi-SNP hotspot
        snps_in_hotspot = set()
        for s in hotspot["snps"]:
            snps_in_hotspot.update(str(s).strip().split("|"))

        # All SNPs not in any hotspot get a single-SNP row so we don't lose them
        all_snp_ids = set(snp_df["snp_id"].astype(str))
        snps_not_in_hotspot = sorted(all_snp_ids - snps_in_hotspot)

        single_rows = []
        for snp_id in snps_not_in_hotspot:
            row = snp_df[snp_df["snp_id"] == snp_id].iloc[0]
            chr_name = row["chr"]
            pos = int(row["pos"])
            single_rows.append({
                "cluster_id": f"{chr_name}_{pos}_{pos}",
                "chr": chr_name,
                "start_bp": pos,
                "end_bp": pos,
                "n_snps": 1,
                "mean_jaccard": "",
                "mean_shared_trans_genes": "",
                "overlap_pvalue": "",
                "snps": snp_id,
            })

        if single_rows:
            single_df = pd.DataFrame(single_rows)
            combined = pd.concat([hotspot, single_df], ignore_index=True)
        else:
            combined = hotspot

        out_name = f"hotspot_groups_winter_{distance_kb}kb.csv"
        combined.to_csv(os.path.join(RESULTS_DIR, out_name), sep=",", index=False)
        print(f"  Hotspot groups (overlap p < {MAX_OVERLAP_PVALUE}): {len(hotspot)}")
        print(f"  Single-SNP rows added (not in any hotspot): {len(single_rows)}")
        print(f"  Total rows in table: {len(combined)}")

    print("\nDone. Results written to", RESULTS_DIR)


if __name__ == "__main__":
    main()
