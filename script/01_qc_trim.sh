#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# Default Parameters (Can be overridden by environment variables)
SAMPLE_TSV="${SAMPLE_TSV:-./sample.tsv}"
RESULTS_DIR="${RESULTS_DIR:-./results}"
QC_DIR="${RESULTS_DIR}/00_qc"
TRIM_DIR="${RESULTS_DIR}/01_trimmed"

# Parallel resource parameters
THREADS="${THREADS:-8}"
FASTQC_THREADS="${FASTQC_THREADS:-${THREADS}}"
TRIM_CONCURRENCY="${TRIM_CONCURRENCY:-8}"
TRIM_THREADS_PER_JOB="${TRIM_THREADS_PER_JOB:-8}"

# V3-V4 Primers (Standard Illumina 16S V3-V4)
PRIMER_F="${PRIMER_F:-CCTACGGGNGGCWGCAG}"
PRIMER_R="${PRIMER_R:-GACTACHVGGGTATCTAATCC}"

echo "========================================================================="
echo "        16S rRNA Quality Control & Primer Trimming (Parallel)            "
echo "========================================================================="
echo "Sample TSV:             ${SAMPLE_TSV}"
echo "Results Directory:      ${RESULTS_DIR}"
echo "FastQC Threads:         ${FASTQC_THREADS}"
echo "Cutadapt Concurrency:   ${TRIM_CONCURRENCY} concurrent jobs"
echo "Threads Per Cutadapt:   ${TRIM_THREADS_PER_JOB} threads/job"
echo "Total Trim Cores:       $((TRIM_CONCURRENCY * TRIM_THREADS_PER_JOB))"
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
raw_fastqs=$(tail -n +2 "${SAMPLE_TSV}" | cut -f2,3 | tr '\t' '\n' | sort -u)

if [ -z "${raw_fastqs}" ]; then
    echo "ERROR: No raw FASTQ paths found in ${SAMPLE_TSV}."
    exit 1
fi

fastqc -t "${FASTQC_THREADS}" -o "${QC_DIR}/raw" ${raw_fastqs}

# Step 2: Primer Trimming with Cutadapt
echo ">>> Step 2: Trimming V3-V4 primers using Cutadapt..."

# Export variables so that they are accessible inside GNU Parallel/xargs subshells
export PRIMER_F PRIMER_R TRIM_DIR TRIM_THREADS_PER_JOB

trim_sample() {
    local sample_name="$1"
    local forward_path="$2"
    local reverse_path="$3"
    
    [ -z "${sample_name}" ] && return 0
    
    if [ ! -f "${forward_path}" ] || [ ! -f "${reverse_path}" ]; then
        echo "Warning: Files for ${sample_name} not found, skipping..."
        return 0
    fi
    
    echo "Processing Sample in Parallel: ${sample_name}"
    
    out_r1="${TRIM_DIR}/${sample_name}_trimmed_R1.fastq.gz"
    out_r2="${TRIM_DIR}/${sample_name}_trimmed_R2.fastq.gz"
    
    if [ -f "${out_r1}" ] && [ -f "${out_r2}" ] && [ -s "${out_r1}" ] && [ -s "${out_r2}" ]; then
        echo "⏭️ Sample ${sample_name} already trimmed and verified. Skipping..."
        return 0
    fi
    
    cutadapt \
        -g "${PRIMER_F}" \
        -G "${PRIMER_R}" \
        -o "${out_r1}" \
        -p "${out_r2}" \
        --discard-untrimmed \
        --minimum-length 200 \
        --cores "${TRIM_THREADS_PER_JOB}" \
        "${forward_path}" \
        "${reverse_path}" \
        > "${TRIM_DIR}/${sample_name}_cutadapt.log" 2>&1
}
export -f trim_sample

# Execute primer trimming in parallel using GNU Parallel (or xargs as a robust fallback)
if command -v parallel &> /dev/null; then
    echo "Running with GNU Parallel..."
    tail -n +2 "${SAMPLE_TSV}" | cut -f1,2,3 | tr '\t' '\n' | parallel -N 3 -j "${TRIM_CONCURRENCY}" trim_sample {1} {2} {3}
else
    echo "Running with xargs fallback..."
    tail -n +2 "${SAMPLE_TSV}" | cut -f1,2,3 | tr '\t' '\n' | xargs -P "${TRIM_CONCURRENCY}" -n 3 bash -c 'trim_sample "$0" "$1" "$2"'
fi

# Step 3: Quality Control of Trimmed Reads
echo ">>> Step 3: Running FastQC on trimmed reads..."
trimmed_fastqs=$(find "${TRIM_DIR}" -maxdepth 1 -name "*.fastq.gz" -o -name "*.fq.gz" | sort)
fastqc -t "${FASTQC_THREADS}" -o "${QC_DIR}/trimmed" ${trimmed_fastqs}

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
