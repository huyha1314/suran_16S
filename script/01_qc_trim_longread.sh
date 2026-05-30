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

# Filtlong Parameters (Length and Quality Thresholds)
FILTLONG_MIN_LEN="${FILTLONG_MIN_LEN:-1000}"
FILTLONG_MIN_Q="${FILTLONG_MIN_Q:-10}"

# Full-length 16S Primers (e.g. 27F and 1492R)
PRIMER_F_LONG="${PRIMER_F_LONG:-AGAGTTTGATCMTGGCTCAG}" # 27F
PRIMER_R_LONG="${PRIMER_R_LONG:-CGGTTACCTTGTTACGACTT}" # 1492R

echo "========================================================================="
echo "        16S Long-Read Quality Control & Primer Trimming (Parallel)       "
echo "========================================================================="
echo "Sample TSV:             ${SAMPLE_TSV}"
echo "Results Directory:      ${RESULTS_DIR}"
echo "FastQC/NanoPlot Cores:  ${FASTQC_THREADS}"
echo "Cutadapt Concurrency:   ${TRIM_CONCURRENCY} concurrent jobs"
echo "Threads Per Job:        ${TRIM_THREADS_PER_JOB} threads/job"
echo "Filtlong Min Length:    ${FILTLONG_MIN_LEN} bp"
echo "Filtlong Min Quality:   Q${FILTLONG_MIN_Q}"
echo "Long-Read Forward:      ${PRIMER_F_LONG}"
echo "Long-Read Reverse:      ${PRIMER_R_LONG}"
echo "========================================================================="

# Create necessary directories
mkdir -p "${QC_DIR}/raw"
mkdir -p "${QC_DIR}/trimmed"
mkdir -p "${QC_DIR}/nanoplot_raw"
mkdir -p "${QC_DIR}/nanoplot_trimmed"
mkdir -p "${TRIM_DIR}"

# Validate sample.tsv exists
if [ ! -f "${SAMPLE_TSV}" ]; then
    echo "ERROR: Standardized sample sheet not found at: ${SAMPLE_TSV}"
    echo "Please generate it first."
    exit 1
fi

# Step 1: Quality Control of Raw Reads
echo ">>> Step 1.1: Running FastQC on raw long reads..."
raw_fastqs=$(tail -n +2 "${SAMPLE_TSV}" | cut -f2 | sort -u)

if [ -z "${raw_fastqs}" ]; then
    echo "ERROR: No raw FASTQ paths found in ${SAMPLE_TSV}."
    exit 1
fi

fastqc -t "${FASTQC_THREADS}" -o "${QC_DIR}/raw" ${raw_fastqs}

echo ">>> Step 1.2: Running NanoPlot on raw long reads..."
NanoPlot -t "${FASTQC_THREADS}" --fastq ${raw_fastqs} -o "${QC_DIR}/nanoplot_raw"

# Step 2: Quality Filtering & Primer Trimming
echo ">>> Step 2: Filtering with Filtlong & Trimming primers using Cutadapt..."

export FILTLONG_MIN_LEN FILTLONG_MIN_Q PRIMER_F_LONG PRIMER_R_LONG TRIM_DIR TRIM_THREADS_PER_JOB

trim_sample_long() {
    local sample_name="$1"
    local forward_path="$2"
    
    [ -z "${sample_name}" ] && return 0
    [ -z "${forward_path}" ] && return 0
    
    if [ ! -f "${forward_path}" ]; then
        echo "Warning: File for ${sample_name} not found at ${forward_path}, skipping..."
        return 0
    fi
    
    echo "Processing Long-Read Sample in Parallel: ${sample_name}"
    
    out_fq="${TRIM_DIR}/${sample_name}_trimmed.fastq.gz"
    
    if [ -f "${out_fq}" ] && [ -s "${out_fq}" ]; then
        echo "⏭️ Sample ${sample_name} already filtered and trimmed. Skipping..."
        return 0
    fi
    
    # 2.1 Quality and Length filtering (Deduplicating on-the-fly to prevent Filtlong duplicate read errors)
    dedup_fq="${TRIM_DIR}/${sample_name}_dedup.fastq"
    gunzip -c "${forward_path}" | awk '
      NR % 4 == 1 { header = $1; line1 = $0; next }
      NR % 4 == 2 { line2 = $0; next }
      NR % 4 == 3 { line3 = $0; next }
      NR % 4 == 0 {
        if (!seen[header]++) {
          print line1
          print line2
          print line3
          print $0
        }
      }
    ' > "${dedup_fq}"

    filt_fq="${TRIM_DIR}/${sample_name}_filt.fastq"
    filtlong \
        --min_length "${FILTLONG_MIN_LEN}" \
        --min_mean_q "${FILTLONG_MIN_Q}" \
        "${dedup_fq}" \
        > "${filt_fq}" 2> "${TRIM_DIR}/${sample_name}_filtlong.log"
        
    # Clean up temp file
    rm -f "${dedup_fq}"
        
    # 2.2 Primer trimming with Cutadapt
    cutadapt \
        -g "${PRIMER_F_LONG}" \
        -a "${PRIMER_R_LONG}" \
        --revcomp \
        -o "${out_fq}" \
        --minimum-length "${FILTLONG_MIN_LEN}" \
        --cores "${TRIM_THREADS_PER_JOB}" \
        "${filt_fq}" \
        > "${TRIM_DIR}/${sample_name}_cutadapt.log" 2>&1
        
    # Clean up temporary uncompressed filt FASTQ file
    rm -f "${filt_fq}"
}
export -f trim_sample_long

# Execute primer trimming in parallel
if command -v parallel &> /dev/null; then
    echo "Running with GNU Parallel..."
    tail -n +2 "${SAMPLE_TSV}" | cut -f1,2 | tr '\t' '\n' | parallel -N 2 -j "${TRIM_CONCURRENCY}" trim_sample_long {1} {2}
else
    echo "Running with xargs fallback..."
    tail -n +2 "${SAMPLE_TSV}" | cut -f1,2 | tr '\t' '\n' | xargs -P "${TRIM_CONCURRENCY}" -n 2 bash -c 'trim_sample_long "$0" "$1"'
fi

# Step 3: Quality Control of Trimmed/Filtered Reads
echo ">>> Step 3.1: Running FastQC on trimmed/filtered long reads..."
trimmed_fastqs=$(find "${TRIM_DIR}" -maxdepth 1 -name "*_trimmed.fastq.gz" | sort)
if [ -n "${trimmed_fastqs}" ]; then
    fastqc -t "${FASTQC_THREADS}" -o "${QC_DIR}/trimmed" ${trimmed_fastqs}
fi

echo ">>> Step 3.2: Running NanoPlot on final trimmed/filtered long reads..."
if [ -n "${trimmed_fastqs}" ]; then
    NanoPlot -t "${FASTQC_THREADS}" --fastq ${trimmed_fastqs} -o "${QC_DIR}/nanoplot_trimmed"
fi

# Step 4: MultiQC Summary Report
echo ">>> Step 4: Generating MultiQC report..."
multiqc -o "${RESULTS_DIR}" \
        -n multiqc_report.html \
        -f \
        "${QC_DIR}/raw" \
        "${TRIM_DIR}" \
        "${QC_DIR}/trimmed"

echo "========================================================================="
echo "Long-Read Quality Control, Filtering and Trimming Completed Successfully!"
echo "Trimmed FASTQ files are saved in:    ${TRIM_DIR}"
echo "FastQC/NanoPlot reports saved in:    ${QC_DIR}"
echo "MultiQC summary is saved in:         ${RESULTS_DIR}/multiqc_report.html"
echo "========================================================================="
