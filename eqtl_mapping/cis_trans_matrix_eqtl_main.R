# cis_trans_matrix_eqtl_main.R

library(MatrixEQTL)
useModel <- modelLINEAR

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  arg_list <- list()
  for (i in seq(1, length(args), by = 2)) {
    key <- gsub("^--", "", args[i])
    value <- args[i + 1]
    arg_list[[key]] <- value
  }
  return(arg_list)
}

params <- parse_args(args)

genotype_file <- params$genotype_file
expression_file <- params$expression_file


#covariate_file <- params$covariate_file
#output_dir <- params$output_dir

dir_path <- params$dir_path
snp_location_file_name <- params$snp_location_file_name
genepos_file_name <- params$genepos_file_name
pvOutputThreshold_cis <- as.numeric(params$pvOutputThreshold_cis)
pvOutputThreshold_trans <- as.numeric(params$pvOutputThreshold_trans)



#errorCovariance <- 1
#errorCovariance <- as.numeric(params$errorCovariance)
#errorCovariance <- params$errorCovariance
cisDist <- params$cisDist

#dir_path <- file.path(dir_path, paste0(dir_name))#"/mnt/winterprojectceph/eqtl_runs/results/run1"
output_file_name_cis <- file.path(dir_path, paste0("cis_eqtl_results_",".txt"))
output_file_name_trans <- file.path(dir_path, paste0("trans_eqtl_results_",".txt"))
log_file <- file.path(dir_path, paste0("log_file_eqtl_",".txt"))
hist_file <- file.path(dir_path, paste0("eqtl_pvalue_histogram_",".png"))

cat("expression_file path:", expression_file, "\n")
cat("genotype_file path:", genotype_file, "\n")

# Load genotype data
snps <- SlicedData$new()
snps$fileDelimiter <- "\t"
snps$fileOmitCharacters <- "NA"
snps$fileSkipRows <- 1
snps$fileSkipColumns <- 1
snps$fileSliceSize <- 2000
snps$LoadFile(genotype_file)


# Load gene expression data
gene <- SlicedData$new()
gene$fileDelimiter <- "\t"
gene$fileOmitCharacters <- "NA"
gene$fileSkipRows <- 1
gene$fileSkipColumns <- 1
gene$fileSliceSize <- 2000
gene$LoadFile(expression_file)


# # Load covariates (empty)
cvrt <- SlicedData$new()
# cvrt$fileDelimiter <- "\t"
# cvrt$fileOmitCharacters <- "NA"
# cvrt$fileSkipRows <- 1
# cvrt$fileSkipColumns <- 1
# cvrt$fileSliceSize <- 2000
# cvrt$LoadFile(covariate_file)


# Load SNP and gene positions
snpspos <- read.table(snp_location_file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names=1)
genepos <- read.table(genepos_file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names=1)


snpspos$SNP <- rownames(snpspos)
snpspos <- snpspos[, c("SNP", "CHR", "POS")]
genepos$geneid <- rownames(genepos)
genepos <- genepos[, c("geneid", "chr", "left", "right")]

snpspos$POS <- as.numeric(snpspos$POS)
genepos$left <- as.numeric(genepos$left)
genepos$right <- as.numeric(genepos$right)


snpspos$CHR <- as.character(snpspos$CHR)
genepos$chr <- as.character(genepos$chr)

# cvrt$columnNames <- gsub("-", ".", cvrt$columnNames)

snps_ids <- snps$columnNames
gene_ids <- gene$columnNames
# cvrt_ids <- cvrt$columnNames
snps$columnNames <- gsub("\"", "", snps$columnNames)
gene$columnNames <- gsub("\"", "", gene$columnNames)

# cvrt$columnNames <- gsub("\"", "", cvrt$columnNames)

print("testing matrixeqtl")

# Run Matrix eQTL
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_trans,
  pvOutputThreshold = pvOutputThreshold_trans,
  useModel = useModel,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
  #noFDRsaveMemory = TRUE
)
# Save cis and trans results
write.table(me$cis$eqtls, file = output_file_name_cis, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(me$trans$eqtls, file = output_file_name_trans, sep = "\t", quote = FALSE, row.names = FALSE)

# Log results
sink(log_file)
cat('Analysis done in:', me$time.in.sec, 'seconds\n')
cat('Detected eQTLs:\n')
#print(me$all$eqtls)
sink()

# Plot histogram
png(hist_file, width = 800, height = 600)
plot(me)
dev.off()

