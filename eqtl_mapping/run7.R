result_dir <- "/mnt/winterprojectceph/eqtl_runs/results"
args <- commandArgs(trailingOnly = FALSE)

#Look for the --file= argument which contains the script path
script_path <- sub("--file=", "", args[grep("--file=", args)])

# Extract just the file name
script_name <- basename(script_path)

print(script_name)

#Create a directory named after the script (without extension)
dir_name <- sub("\\.R$", "", script_name)  # Remove .R extension
dir_path <- file.path(result_dir, dir_name)   # Create full path in current working directory

#Use system command to create the directory
system(paste("mkdir -p", shQuote(dir_path)))

#Confirm creation
cat("Directory created at:", dir_path, "\n")



#run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")  # e.g., "20250918_1206"
#cat("Run ID:", run_id, "\n")
run_id = ""
#Input files Path
genotype_file <- "/mnt/winterprojectceph//vcf_file/genotype_matrix_output1.txt"
expression_file <- "/mnt/winterprojectceph/winter_pipe/combine_batch_results/gf_cnt_mtx/results/gene_count_matrix_summer_vsd_renamed.csv"
covariate_file <- "/mnt/winterprojectceph/vcf_file/covariates_for_matrix_eqtl.txt"

genotype_aligned <- file.path(dir_path, "genotype_matrix_new_summer_aligned.txt")
expression_aligned <- file.path(dir_path, "expression_matrix_newsummeraligned.txt")
covariates_aligned <- file.path(dir_path, "covariates_matrix_new_summeraligned.txt")

#snp_location_file_name <- "/mnt/winterprojectceph/vcf_file/plink_new/snp_location_output_1000.txt"
snp_location_file_name <- "/mnt/winterprojectceph/vcf_file/snp_location_test1.txt"
genepos_file_name <- "/mnt/winterprojectceph/vcf_file/genepos_clean.txt"

#Output files Path

output_file_name_cis <- file.path(dir_path, paste0("cis_eqtl_results_full_1e-5",".txt"))
output_file_name_trans <- file.path(dir_path, paste0("trans_eqtl_results_1e-5",".txt"))
log_file <- file.path(dir_path, paste0("log_file_eqtl_1e-5",".txt"))
hist_file <- file.path(dir_path, paste0("eqtl_pvalue_histogram_1e-5",".png"))


#Parameters

pvOutputThreshold_cis <- 1
pvOutputThreshold_trans <- 1
#errorCovariance <- numeric()

cisDist <- 1e6


 

system2("Rscript", args = c("align_matrix.R",
                            "--genotype_file", genotype_file,
                            "--expression_file", expression_file,
			    			"--covariate_file", covariate_file,
			    			"--output_dir", dir_path,
			    			"--run_id", run_id))

# system2("Rscript", args = c("cis_trans_matrix_eqtl_main.R",
#                             "--genotype_file", genotype_aligned,
#                             "--expression_file", expression_aligned,
# 			    			"--covariate_file", covariates_aligned,
# 			    			"--dir_path", dir_path,
# 			    			"--snp_location_file_name", snp_location_file_name,
# 							"--genepos_file_name", genepos_file_name,
# 							"--pvOutputThreshold_cis", pvOutputThreshold_cis,
# 							"--pvOutputThreshold_trans", pvOutputThreshold_trans,
# 							"--cisDist", cisDist
# 							))


#"--errorCovariance", errorCovariance,
