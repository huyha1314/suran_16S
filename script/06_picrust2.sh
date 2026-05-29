#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# Parameters
DADA2_DIR="${DADA2_DIR:-./results/02_dada2}"
PHYLO_DIR="${PHYLO_DIR:-./results/03_phyloseq}"
PICRUSt2_DIR="${PICRUSt2_DIR:-./results/06_picrust2}"
THREADS="${THREADS:-8}"

# Input files
FASTA_IN="${DADA2_DIR}/asvs.fasta"
TABLE_IN="${PHYLO_DIR}/asv_table_picrust.tsv"

echo "========================================================================="
echo "               16S rRNA PICRUSt2 Functional Prediction                   "
echo "========================================================================="
echo "ASV Fasta:       ${FASTA_IN}"
echo "ASV Count Table: ${TABLE_IN}"
echo "Output Directory: ${PICRUSt2_DIR}"
echo "Threads:         ${THREADS}"
echo "========================================================================="

# Check if input files exist
if [ ! -f "${FASTA_IN}" ] || [ ! -f "${TABLE_IN}" ]; then
    echo "ERROR: Required input files not found. Ensure DADA2 and Phyloseq steps are completed."
    echo "Fasta exists: $([ -f "${FASTA_IN}" ] && echo "Yes" || echo "No")"
    echo "Table exists: $([ -f "${TABLE_IN}" ] && echo "Yes" || echo "No")"
    exit 1
fi

# Clean up existing output directory to prevent PICRUSt2 from stopping
rm -rf "${PICRUSt2_DIR}"

# Convert to absolute paths so PICRUSt2 R subprocesses can find files
# regardless of their working directory
ABS_FASTA_IN="$(realpath "${FASTA_IN}")"
ABS_TABLE_IN="$(realpath "${TABLE_IN}")"
ABS_PICRUST2_DIR="$(cd "$(dirname "${PICRUSt2_DIR}")" && pwd)/$(basename "${PICRUSt2_DIR}")"

# Run PICRUSt2 pipeline
echo ">>> Running PICRUSt2 pipeline..."
picrust2_pipeline.py \
    -s "${ABS_FASTA_IN}" \
    -i "${ABS_TABLE_IN}" \
    -o "${ABS_PICRUST2_DIR}" \
    -p "${THREADS}"

# Decompress primary output tables for direct user visibility and downstream R analyses
echo ">>> Decompressing predicted functional abundance tables..."

if [ -f "${ABS_PICRUST2_DIR}/pathways_out/path_abun_unstrat.tsv.gz" ]; then
    gunzip -f "${ABS_PICRUST2_DIR}/pathways_out/path_abun_unstrat.tsv.gz"
    echo "Unstratified Pathway Abundance: ${ABS_PICRUST2_DIR}/pathways_out/path_abun_unstrat.tsv"
fi

if [ -f "${ABS_PICRUST2_DIR}/ec_metagenome_out/pred_metagenome_unstrat.tsv.gz" ]; then
    gunzip -f "${ABS_PICRUST2_DIR}/ec_metagenome_out/pred_metagenome_unstrat.tsv.gz"
    echo "Unstratified EC Metagenome:    ${ABS_PICRUST2_DIR}/ec_metagenome_out/pred_metagenome_unstrat.tsv"
fi

if [ -f "${ABS_PICRUST2_DIR}/ko_metagenome_out/pred_metagenome_unstrat.tsv.gz" ]; then
    gunzip -f "${ABS_PICRUST2_DIR}/ko_metagenome_out/pred_metagenome_unstrat.tsv.gz"
    echo "Unstratified KO Metagenome:    ${ABS_PICRUST2_DIR}/ko_metagenome_out/pred_metagenome_unstrat.tsv"
fi

echo "========================================================================="
echo "Functional prediction completed successfully!"
echo "Outputs saved in: ${ABS_PICRUST2_DIR}"
echo "========================================================================="
