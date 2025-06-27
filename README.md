# WINTR
This is the pipeline to process RNA sequencing data.

**main.nf**
This Nextflow workflow performs a complete RNA-seq data processing and quantification pipeline, beginning with raw interleaved FASTQ files and producing gene- and transcript-level TPM and count matrices. The pipeline starts by splitting interleaved reads, followed by quality control using FastQC before and after trimming with Trimmomatic. Trimmed paired-end reads are aligned to a reference genome using HISAT2, and the resulting BAM files are used for transcript assembly and quantification with StringTie. Each sample’s StringTie output (GTF and gene expression .ctab file) is further processed to extract TPM values and merged into comprehensive matrices. The pipeline includes custom Python scripts to annotate TPM values in .ctab files and compile sample lists used by prepDE.py to generate count matrices. Finally, it creates a gene- and transcript-level TPM matrix, providing outputs suitable for downstream expression analysis. 

**nextflow.config**
It specifies file paths for input data (like FASTQ files, reference genome, annotation files) and tools (scripts and HISAT2 index). Parameters like number of threads and script locations are centralized here for consistency.



**add_tpm_to_ctab.py**
This script extracts TPM values for each transcript from a GTF file and adds them to the corresponding t_data.ctab file based on matching transcript_id (t_name). It reads both files, merges them on the t_name column, and creates a new output file with the added TPM column. The merging ensures each transcript entry in ctab gets the correct TPM from the GTF. Finally, it saves the merged file as a tab-separated output.


**build_tpm_matrices.py**
This Python script generates transcript- and gene-level TPM expression matrices from a list of .ctab files. It reads each file, verifies required columns (TPM, t_name, and gene_id), and aggregates data by t_name and gene_id. The TPM values are combined across all samples and saved as two .tsv files. Errors and warnings are printed for missing data or formatting issues.


**prepDE.py**

This script, originally from the StringTie repository (https://github.com/gpertea/stringtie/blob/master/prepDE.py), generates gene- and transcript-level read count matrices from a set of GTF files produced by stringtie -e. It parses coverage data and calculates counts based on transcript lengths and read coverage, supporting optional clustering of overlapping genes. The results are written to two CSV files (gene_count_matrix.csv and transcript_count_matrix.csv). 

