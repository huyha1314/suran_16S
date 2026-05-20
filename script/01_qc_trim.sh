#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# Default Parameters (Can be overridden by environment variables)
SAMPLE_TSV="${SAMPLE_TSV:-./data/sample.tsv}"
RESULTS_DIR="${RESULTS_DIR:-./results}"
QC_DIR="${RESULTS_DIR}/00_qc"
TRIM_DIR="${RESULTS_DIR}/01_trimmed"
THREADS="${THREADS:-8}"

# V3-V4 Primers (Standard Illumina 16S V3-V4)
PRIMER_F="${PRIMER_F:-CCTACGGGNGGCWGCAG}"
PRIMER_R="${PRIMER_R:-GACTACHVGGGTATCTAATCC}"

echo "========================================================================="
echo "               16S rRNA V3-V4 Quality Control & Primer Trimming          "
echo "========================================================================="
echo "Sample TSV:       ${SAMPLE_TSV}"
echo "Results Directory: ${RESULTS_DIR}"
echo "Threads:          ${THREADS}"
echo "Forward Primer:   ${PRIMER_F}"
echo "Reverse Primer:   ${PRIMER_R}"
echo "========================================================================="

# Create necessary directories
mkdir -p "${QC_DIR}/raw"
mkdir -p "${QC_DIR}/trimmed"
mkdir -p "${TRIM_DIR}"

# Validate sample.tsv exists
if [ ! -f "${SAMPLE_TSV}" ]; then
    echo "ERROR: Standardized sample sheet not found at: ${SAMPLE_TSV}"
    echo "Please generate it first by running: pixi run generate_samples"
    exit 1
fi

# Step 1: Quality Control of Raw Reads
echo ">>> Step 1: Running FastQC on raw reads..."
# Extract unique fastq paths from sample sheet columns 2 (ForwardPath) and 3 (ReversePath)
raw_fastqs=$(tail -n +2 "${SAMPLE_TSV}" | cut -f2,3 | tr '\t' '\n' | sort -u)

if [ -z "${raw_fastqs}" ]; then
    echo "ERROR: No raw FASTQ paths found in ${SAMPLE_TSV}."
    exit 1
fi

fastqc -t "${THREADS}" -o "${QC_DIR}/raw" ${raw_fastqs}

# Step 2: Primer Trimming with Cutadapt
echo ">>> Step 2: Trimming V3-V4 primers using Cutadapt..."

tail -n +2 "${SAMPLE_TSV}" | while IFS=$'\t' read -r sample_name forward_path reverse_path; do
    # Skip empty lines
    [ -z "${sample_name}" ] && continue
    
    echo "Processing Sample: ${sample_name}"
    
    if [ ! -f "${forward_path}" ] || [ ! -f "${reverse_path}" ]; then
        echo "Warning: Files for ${sample_name} not found, skipping..."
        continue
    fi
    
    # Define outputs using trimmed nomenclature
    out_r1="${TRIM_DIR}/${sample_name}_trimmed_R1.fastq.gz"
    out_r2="${TRIM_DIR}/${sample_name}_trimmed_R2.fastq.gz"
    
    # Run Cutadapt
    cutadapt \
        -g "${PRIMER_F}" \
        -G "${PRIMER_R}" \
        -o "${out_r1}" \
        -p "${out_r2}" \
        --discard-untrimmed \
        --minimum-length 200 \
        --cores "${THREADS}" \
        "${forward_path}" \
        "${reverse_path}" \
        > "${TRIM_DIR}/${sample_name}_cutadapt.log"
done

# Step 3: Quality Control of Trimmed Reads
echo ">>> Step 3: Running FastQC on trimmed reads..."
trimmed_fastqs=$(find "${TRIM_DIR}" -maxdepth 1 -name "*.fastq.gz" -o -name "*.fq.gz" | sort)
fastqc -t "${THREADS}" -o "${QC_DIR}/trimmed" ${trimmed_fastqs}

# Step 4: MultiQC Summary Report
echo ">>> Step 4: Generating MultiQC report..."
multiqc -o "${RESULTS_DIR}" \
        -n multiqc_report.html \
        -f \
        "${QC_DIR}/raw" \
        "${TRIM_DIR}" \
        "${QC_DIR}/trimmed"

echo "========================================================================="
echo "Quality Control and Primer Trimming Completed Successfully!"
echo "Trimmed FASTQ files are saved in: ${TRIM_DIR}"
echo "QC reports are saved in:          ${QC_DIR}"
echo "MultiQC summary is saved in:      ${RESULTS_DIR}/multiqc_report.html"
echo "========================================================================="
