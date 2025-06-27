import pandas as pd
import argparse

def build_tpm_matrices_from_sample_list(sample_list_file, out_prefix):
    transcript_dfs = []
    gene_dfs = []

    with open(sample_list_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            try:
                sample_name, file_path = line.strip().split('\t')
            except ValueError:
                print(f"[ERROR] Malformed line: {line.strip()}")
                continue

            print(f"[INFO] Processing {sample_name} → {file_path}")
            try:
                df = pd.read_csv(file_path, sep='\t')
            except Exception as e:
                print(f"[ERROR] Could not read {file_path}: {e}")
                continue

            df.columns = df.columns.str.strip()

            missing_cols = [c for c in ['TPM', 't_name', 'gene_id'] if c not in df.columns]
            if missing_cols:
                print(f"[WARNING] Skipping {sample_name}: missing columns {missing_cols}")
                continue

            # ✅ Transcript-level TPM — use t_name
            tx_df = df[['t_name', 'TPM']].copy()
            tx_df.rename(columns={'TPM': sample_name}, inplace=True)
            tx_df.set_index('t_name', inplace=True)
            tx_df.fillna(0, inplace=True)
            transcript_dfs.append(tx_df)

            # ✅ Gene-level TPM — use gene_id
            gene_df = df[['gene_id', 'TPM']].copy()
            gene_df = gene_df.groupby('gene_id', as_index=False).sum()
            gene_df.rename(columns={'TPM': sample_name}, inplace=True)
            gene_df.set_index('gene_id', inplace=True)
            gene_df.fillna(0, inplace=True)
            gene_dfs.append(gene_df)

    if transcript_dfs:
        transcript_matrix = pd.concat(transcript_dfs, axis=1)
        transcript_matrix.fillna(0, inplace=True)
        transcript_matrix.to_csv(f'{out_prefix}_transcript_matrix.tsv', sep='\t')
        print(f"[✓] Transcript-level TPM matrix written to: {out_prefix}_transcript_matrix.tsv")
    else:
        print("[ERROR] No transcript-level TPM data was parsed.")

    if gene_dfs:
        gene_matrix = pd.concat(gene_dfs, axis=1)
        gene_matrix.fillna(0, inplace=True)
        gene_matrix.to_csv(f'{out_prefix}_gene_matrix.tsv', sep='\t')
        print(f"[✓] Gene-level TPM matrix written to: {out_prefix}_gene_matrix.tsv")
    else:
        print("[ERROR] No gene-level TPM data was parsed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_list", required=True, help="Path to ctab_sample_lst.txt")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output matrix files")
    args = parser.parse_args()

    build_tpm_matrices_from_sample_list(args.sample_list, args.out_prefix)

