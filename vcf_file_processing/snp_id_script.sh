bcftools annotate \
  -x ID \
  -I +'%CHROM\_%POS' \
  -o fixed_snps.vcf.gz \
  -O z \
  populus_variation_miss010_maf005.vcf.gz

