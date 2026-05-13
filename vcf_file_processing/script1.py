import pandas as pd

# Load the .traw file
traw_file = "populus_output1_genotypes.traw"  # Replace with your actual file path
df = pd.read_csv(traw_file, delim_whitespace=True)

# Extract SNP location info
snp_location = df[['SNP', 'CHR', 'POS']]
snp_location.to_csv("snp_location_output1.txt", sep='\t', index=False)

# Clean sample column names
sample_cols = df.columns[6:]  # First 6 columns are metadata
clean_sample_cols = [col.split('_')[0] for col in sample_cols]

# Extract genotype matrix
genotype_matrix = df[sample_cols].copy()
genotype_matrix.columns = clean_sample_cols
genotype_matrix.insert(0, 'SNP', df['SNP'])  # Add SNP as first column
genotype_matrix.to_csv("genotype_matrix_output1.txt", sep='\t', index=False)

print("Files generated:")
print("- snp_location_output1.txt")
print("- genotype_matrix_output1.txt")
