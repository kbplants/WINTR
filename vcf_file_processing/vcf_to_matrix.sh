#!/bin/bash
set -euo pipefail

# Input gzipped VCF file
VCF_FILE="filtered_final.recode.vcf.gz"
OUT_PREFIX="populus_output1"

# Step 1: Convert gzipped VCF to PLINK binary format
plink --vcf "$VCF_FILE" \
      --make-bed \
      --allow-extra-chr \
      --out "$OUT_PREFIX"

# Step 2: Replace missing SNP IDs with unique ones (if BIM contains '.')
awk 'BEGIN{OFS="\t"} {if($2==".") $2="snp"NR; print}' "${OUT_PREFIX}.bim" > "${OUT_PREFIX}.bim.tmp"
mv "${OUT_PREFIX}.bim.tmp" "${OUT_PREFIX}.bim"

# Step 3: Extract genotype matrix in raw format (transpose = samples as rows, SNPs as columns)
plink --bfile "$OUT_PREFIX" \
      --recode A-transpose \
      --allow-extra-chr \
      --out "${OUT_PREFIX}_genotypes"

# Step 4: Extract SNP location file from updated BIM
awk 'BEGIN{OFS="\t"} {print $2, $1, $4}' "${OUT_PREFIX}.bim" > "${OUT_PREFIX}_snp_locations.txt"

echo "Done! Outputs:"
echo "  - PLINK binary files: ${OUT_PREFIX}.bed, ${OUT_PREFIX}.bim, ${OUT_PREFIX}.fam"
echo "  - Genotype matrix: ${OUT_PREFIX}_genotypes.raw"
echo "  - SNP location file: ${OUT_PREFIX}_snp_locations.txt"
