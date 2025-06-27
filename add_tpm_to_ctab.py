import pandas as pd
import re
import argparse

def extract_tpm_from_gtf(gtf_path):
    """Extract transcript_id → TPM from GTF."""
    tpm_records = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or "\ttranscript\t" not in line:
                continue
            fields = line.strip().split('\t')
            attr = fields[8]
            tid_match = re.search(r'transcript_id "([^"]+)"', attr)
            tpm_match = re.search(r'TPM "([^"]+)"', attr)
            if tid_match and tpm_match:
                tid = tid_match.group(1)
                tpm = float(tpm_match.group(1))
                tpm_records.append((tid, tpm))
    return pd.DataFrame(tpm_records, columns=["t_name", "TPM"])

def main(gtf_file, ctab_file, output_file):
    tpm_df = extract_tpm_from_gtf(gtf_file)
    ctab_df = pd.read_csv(ctab_file, sep='\t')

    # Ensure t_name columns match type
    tpm_df['t_name'] = tpm_df['t_name'].astype(str)
    ctab_df['t_name'] = ctab_df['t_name'].astype(str)

    # Merge TPM values into ctab using t_name
    merged = pd.merge(ctab_df, tpm_df, on='t_name', how='left')
    merged.to_csv(output_file, sep='\t', index=False)
    print(f"[✓] Output written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', required=True, help='Input GTF file with TPM')
    parser.add_argument('--ctab', required=True, help='Original t_data.ctab file')
    parser.add_argument('--output', required=True, help='Output file with TPM column')
    args = parser.parse_args()

    main(args.gtf, args.ctab, args.output)

