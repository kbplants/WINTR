#!/usr/bin/env python3
"""
Filter winter trans eQTL results:
  - pvalue < 1e-6
  - Keep only SNPs with >= 50 genes
"""
import os
import pandas as pd

WORKDIR = "/mnt/winterprojectceph/new_winter/001_deseq/test/results/without_threshold/deg_fdr_0.05/024_filter_hotspot"
PVAL_THRESH = 1e-6
MIN_GENES_PER_SNP = 50

INPUT_FILE = os.path.join(WORKDIR, "data", "trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new.txt")
OUTPUT_FILE = os.path.join(WORKDIR, "results", "winter_trans_eqtl_filtered_p1e-6_min50genes.txt")
SUMMARY_FILE = os.path.join(WORKDIR, "results", "filter_summary.txt")


def main():
    os.makedirs(os.path.join(WORKDIR, "results"), exist_ok=True)

    df = pd.read_csv(INPUT_FILE, sep="\t")
    n_total = len(df)
    n_snps_total = df["snps"].nunique()

    # Filter by pvalue
    df_f = df[df["pvalue"] < PVAL_THRESH].copy()
    n_after_pval = len(df_f)
    n_snps_after_pval = df_f["snps"].nunique()

    # Count unique genes per SNP
    genes_per_snp = df_f.groupby("snps")["gene"].nunique()
    snps_keep = genes_per_snp[genes_per_snp >= MIN_GENES_PER_SNP].index.tolist()

    # Keep only rows for SNPs with >= 50 genes
    df_out = df_f[df_f["snps"].isin(snps_keep)]
    n_final = len(df_out)
    n_snps_final = df_out["snps"].nunique()

    df_out.to_csv(OUTPUT_FILE, sep="\t", index=False)

    summary_lines = [
        "Winter trans eQTL filter summary",
        "================================",
        f"Input file: {INPUT_FILE}",
        f"Filters: pvalue < {PVAL_THRESH}, SNPs with >= {MIN_GENES_PER_SNP} genes",
        "",
        f"Total associations (input):     {n_total}",
        f"Unique SNPs (input):           {n_snps_total}",
        f"Associations after pvalue:     {n_after_pval}",
        f"Unique SNPs after pvalue:      {n_snps_after_pval}",
        f"Associations after SNP filter: {n_final}",
        f"Unique SNPs remaining:         {n_snps_final}",
        "",
        f"Output: {OUTPUT_FILE}",
    ]
    summary = "\n".join(summary_lines)
    with open(SUMMARY_FILE, "w") as f:
        f.write(summary)

    print(summary)
    return n_snps_final


if __name__ == "__main__":
    main()
