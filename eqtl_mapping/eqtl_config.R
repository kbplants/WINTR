# eqtl_config.R

# Model type
useModel <- modelLINEAR

# File paths

expression_file_name <- "/mnt/winterprojectceph/eqtl_runs/results/run8/expression_matrix_aligned_cov_r.txt"
SNP_file_name <- "/mnt/winterprojectceph/eqtl_runs/results/run8/genotype_matrix_aligned_cov_r.txt"
covariates_file_name <- "/mnt/winterprojectceph/eqtl_runs/results/run8/covariates_matrix_aligned_cov_r.txt"


#snp_location_file_name <- "/mnt/winterprojectceph/vcf_file/snp_location_output1.txt"
snp_location_file_name <- "/mnt/winterprojectceph/vcf_file/snp_location_output1.txt"
genepos_file_name <- "/mnt/winterprojectceph/vcf_file/genepos_clean.txt"  # Add this if needed



#delete this line as soon as winter data is run

results_dir <- "/mnt/winterprojectceph/eqtl_runs/results/run8"
output_file_name_cis <- file.path(results_dir, "cis_eqtl_winter_bark_xylem_gene_filter_pf_all_new.txt")
output_file_name_trans <- file.path(results_dir, "trans_eqtl_results_winter_bark_xylem_gene_filter_pf_all_new.txt")
log_file <- file.path(results_dir, "log_file_eqtl_winter_bark_xylem_gene_filter_all_new.txt")
hist_file <- file.path(results_dir, "eqtl_pvalue_histogram_winter_bark_xylem_gene_filter_SNP_all_new.png")
####

# Thresholds and settings
#pvOutputThreshold <- 1e-5

###
errorCovariance <- numeric()

cisDist <- 1e6  # 1 Mb window for cis-eQTLs
pvOutputThreshold_cis <- 1e-5
pvOutputThreshold_trans <- 1e-5

#Remove this later
# pvOutputThreshold_cis <- 1
# pvOutputThreshold_trans <- 1

###
