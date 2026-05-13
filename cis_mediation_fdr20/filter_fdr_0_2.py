#!/usr/bin/env python3

import os
import sys
from collections import defaultdict
from typing import Dict, List, Set


def filter_file_by_fdr(path: str, threshold: float = 0.2) -> str:
    """
    Read a tab-delimited file, filter rows by FDR <= threshold,
    and write a new file with '_fdr_20' inserted before the extension
    in the same directory.
    """
    dirname, filename = os.path.split(path)
    name, ext = os.path.splitext(filename)
    if not ext:
        # fall back to appending to whole name if there is no extension
        out_name = f"{name}_fdr_20"
    else:
        out_name = f"{name}_fdr_20{ext}"
    out_path = os.path.join(dirname, out_name)

    with open(path, "r") as fin:
        header = fin.readline()
        if not header:
            # empty file, nothing to filter
            return out_path

        cols = header.rstrip("\n").split("\t")
        try:
            fdr_idx = cols.index("FDR")
        except ValueError:
            # If there is no FDR column (e.g. genotype matrices), skip filtering
            print(f"Warning: No 'FDR' column found in {path}; skipping FDR filter.", file=sys.stderr)
            return path

        with open(out_path, "w") as fout:
            fout.write(header)
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= fdr_idx:
                    continue
                try:
                    fdr_val = float(parts[fdr_idx])
                except ValueError:
                    # skip lines where FDR is not a float
                    continue
                if fdr_val <= threshold:
                    fout.write(line)

    return out_path


def filter_snps_by_gene_count(path: str, min_genes: int = 50) -> str:
    """
    From a tab-delimited file with 'snps' and 'gene' columns, retain only
    rows where the SNP has at least `min_genes` associated (trans) genes.

    Writes a new file with '_atleast50genes' appended before the extension
    and returns the new path.
    """
    dirname, filename = os.path.split(path)
    name, ext = os.path.splitext(filename)
    out_name = f"{name}_atleast50genes{ext}"
    out_path = os.path.join(dirname, out_name)

    with open(path, "r") as fin:
        header = fin.readline()
        if not header:
            return out_path

        cols = header.rstrip("\n").split("\t")
        try:
            snp_idx = cols.index("snps")
            gene_idx = cols.index("gene")
        except ValueError:
            raise RuntimeError(f"Required 'snps' or 'gene' columns not found in {path}")

        snp_to_genes: Dict[str, Set[str]] = defaultdict(set)
        lines: List[str] = []

        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(snp_idx, gene_idx):
                continue
            snp = parts[snp_idx]
            gene = parts[gene_idx]
            snp_to_genes[snp].add(gene)
            lines.append(line)

    # Determine which SNPs meet the threshold
    keep_snps = {snp for snp, genes in snp_to_genes.items() if len(genes) >= min_genes}

    with open(out_path, "w") as fout:
        fout.write(header)
        if not keep_snps:
            return out_path
        for line in lines:
            parts = line.rstrip("\n").split("\t")
            snp = parts[snp_idx]
            if snp in keep_snps:
                fout.write(line)

    return out_path


def main(argv: List[str]) -> None:
    """
    If file paths are provided as arguments, filter only those.
    Otherwise, filter all .txt files in the ../data directory.
    """
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")

    if argv:
        input_files = [os.path.abspath(p) for p in argv]
    else:
        input_files = []
        for name in os.listdir(data_dir):
            if not name.endswith(".txt"):
                continue
            # Skip already filtered / derived files when running in default mode
            if "_fdr_20" in name or "atleast50genes" in name:
                continue
            # Skip genotype matrices and other non-eQTL files (no FDR filtering)
            if "genotype_matrix" in name:
                continue
            input_files.append(os.path.join(data_dir, name))

    if not input_files:
        print("No input files found to filter.", file=sys.stderr)
        return

    print(f"Filtering {len(input_files)} file(s) to FDR <= 0.2...")
    fdr_files: List[str] = []
    for path in input_files:
        print(f"  - {path}")
        fdr_files.append(filter_file_by_fdr(path, threshold=0.2))

    # Additional filtering step for trans-eQTL result files:
    # from the FDR-filtered files, retain only SNPs with at least 50
    # associated trans genes (unique 'gene' values per SNP).
    print("Applying SNP >= 50 trans genes filter for trans-eQTL files...")
    for fdr_path in fdr_files:
        fname = os.path.basename(fdr_path)
        if fname.startswith("trans_eqtl_results_") and "_fdr_20" in fname:
            print(f"  - {fdr_path}")
            filter_snps_by_gene_count(fdr_path, min_genes=50)

    print("Done.")


if __name__ == "__main__":
    main(sys.argv[1:])

