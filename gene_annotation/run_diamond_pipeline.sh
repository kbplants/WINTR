#!/usr/bin/env bash
set -euo pipefail

# DIAMOND pipeline: build a protein DIAMOND database from a reference FASTA,
# then run a blastp search of a query protein FASTA against that database.
# Results are written to the project results/ folder.

PROJECT_DIR="/mnt/winterprojectceph/gene_annot"

# =============================================================================
# EDIT HERE: path to the reference protein FASTA used to BUILD the DIAMOND DB.
# This is the database source (e.g. UniProt plant reference proteins).
# =============================================================================
REFERENCE_FASTA="/mnt/winterprojectceph/gene_annot/data/uniprot/plant_reference_proteins.fasta.gz"

# =============================================================================
# EDIT HERE: path to the QUERY protein FASTA to search against the DIAMOND DB.
# =============================================================================
QUERY_FASTA="/mnt/winterprojectceph/gene_annot/data/populus_trichocarpa/Ptrichocarpa_444_v3.1.protein.lignin_genes.fa.gz"

# Derive the DIAMOND database prefix from the reference FASTA so the .dmnd file
# is created next to it with a matching basename (e.g. plant_reference_proteins
# .fasta.gz -> plant_reference_proteins.dmnd).
ref_basename="$(basename "${REFERENCE_FASTA}")"
ref_basename="${ref_basename%.gz}"
ref_basename="${ref_basename%.fasta}"
ref_basename="${ref_basename%.fa}"
ref_basename="${ref_basename%.faa}"
DB_PREFIX="$(dirname "${REFERENCE_FASTA}")/${ref_basename}"

# Same basename stripping for the query (used in the default hits filename).
query_basename="$(basename "${QUERY_FASTA}")"
query_basename="${query_basename%.gz}"
query_basename="${query_basename%.fasta}"
query_basename="${query_basename%.fa}"
query_basename="${query_basename%.faa}"

# Results directory for DIAMOND search output
OUT_DIR="${PROJECT_DIR}/results/diamond_db"

# Output hits file:
# By default OUTFILE is NOT a fixed name. It is derived from the query and
# reference FASTA file stems so different EDIT_HERE paths produce different
# result files and you do not overwrite previous runs by accident.
#
# Override: set OUTFILE before running, e.g.
#   OUTFILE=/path/to/my_hits.tsv bash scripts/run_diamond_pipeline.sh
#
# If you re-run with the same query and reference and still want a new file,
# set APPEND_TIMESTAMP=1 to add a YYYYMMDD_HHMMSS suffix.
OUTFILE="${OUTFILE:-}"
APPEND_TIMESTAMP="${APPEND_TIMESTAMP:-0}"

THREADS="${THREADS:-8}"
EVALUE="${EVALUE:-1e-5}"
FORCE_REBUILD="${FORCE_REBUILD:-0}"

if ! command -v diamond >/dev/null 2>&1; then
    echo "ERROR: diamond is not available on PATH." >&2
    exit 1
fi

if [[ ! -s "${REFERENCE_FASTA}" ]]; then
    echo "ERROR: reference FASTA not found or empty: ${REFERENCE_FASTA}" >&2
    exit 1
fi

if [[ ! -s "${QUERY_FASTA}" ]]; then
    echo "ERROR: query FASTA not found or empty: ${QUERY_FASTA}" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

if [[ -z "${OUTFILE}" ]]; then
    OUTFILE="${OUT_DIR}/${query_basename}_vs_${ref_basename}_diamond_blastp.tsv"
    if [[ "${APPEND_TIMESTAMP}" == "1" ]]; then
        OUTFILE="${OUT_DIR}/${query_basename}_vs_${ref_basename}_$(date +%Y%m%d_%H%M%S)_diamond_blastp.tsv"
    fi
fi

# Step 1: build the DIAMOND database.
# If "${DB_PREFIX}.dmnd" already exists on disk, the build is skipped so we do
# not re-run "diamond makedb" against the reference. Set FORCE_REBUILD=1 to
# force a rebuild even when the database file already exists.
if [[ -s "${DB_PREFIX}.dmnd" && "${FORCE_REBUILD}" != "1" ]]; then
    echo "DIAMOND database already exists, skipping makedb:"
    echo "  database: ${DB_PREFIX}.dmnd"
    echo "  (set FORCE_REBUILD=1 to rebuild it from ${REFERENCE_FASTA})"
else
    echo "Building DIAMOND database:"
    echo "  input:  ${REFERENCE_FASTA}"
    echo "  output: ${DB_PREFIX}.dmnd"
    diamond makedb --in "${REFERENCE_FASTA}" -d "${DB_PREFIX}"
fi

# Step 2: run DIAMOND blastp of the query against the database
echo "Running DIAMOND blastp:"
echo "  query:  ${QUERY_FASTA}"
echo "  db:     ${DB_PREFIX}.dmnd"
echo "  output: ${OUTFILE}"

diamond blastp \
    -d "${DB_PREFIX}" \
    -q "${QUERY_FASTA}" \
    -o "${OUTFILE}" \
    --threads "${THREADS}" \
    --evalue "${EVALUE}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

echo "Done: ${OUTFILE}"
