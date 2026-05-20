#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status, 
# if an uninitialized variable is used, or if a piped command fails.
set -euo pipefail

# =========================================================================
#      LONGITUDINAL 16S WORKFLOW ORCHESTRATOR & CONFIGURATION PANEL
# =========================================================================

# -------------------------------------------------------------------------
# 1. Global Pipeline Configuration (Edit these parameters as needed)
# -------------------------------------------------------------------------
export THREADS="${THREADS:-8}"                  # Number of CPU cores to allocate

# Core Directory & File Configs
export DATA_DIR="${DATA_DIR:-./data}"           # Path containing raw FASTQ files
export SAMPLE_TSV="${SAMPLE_TSV:-${DATA_DIR}/sample.tsv}" # Standardized sample sheet TSV
export METADATA_PATH="${METADATA_PATH:-${DATA_DIR}/metadata.tsv}" # Standardized metadata TSV
export DB_DIR="${DB_DIR:-./data/db}"           # Directory for reference databases

# Outputs Configs
export RESULTS_DIR="${RESULTS_DIR:-./results}"
export QC_DIR="${QC_DIR:-${RESULTS_DIR}/01_trimmed/qc}"
export TRIM_DIR="${TRIM_DIR:-${RESULTS_DIR}/01_trimmed}"
export DADA2_DIR="${DADA2_DIR:-${RESULTS_DIR}/02_dada2}"
export PHYLO_DIR="${PHYLO_DIR:-${RESULTS_DIR}/03_phyloseq}"
export STATS_DIR="${STATS_DIR:-${RESULTS_DIR}/04_stats}"
export NETWORK_DIR="${NETWORK_DIR:-${RESULTS_DIR}/05_network}"
export PICRUSt2_DIR="${PICRUSt2_DIR:-${RESULTS_DIR}/06_picrust2}"

# Primer Settings (Standard V3-V4 primers)
export PRIMER_F="${PRIMER_F:-CCTACGGGNGGCWGCAG}" # Forward Primer (341F)
export PRIMER_R="${PRIMER_R:-GACTACHVGGGTATCTAATCC}" # Reverse Primer (805R)

# DADA2 Denoising Settings (Asymmetric PE300 Quality Rescuer)
export TRUNC_LEN_F="${TRUNC_LEN_F:-270}"        # Truncate forward reads at 270bp
export TRUNC_LEN_R="${TRUNC_LEN_R:-210}"        # Truncate reverse reads at 210bp (erroneous tails dropped)
export MAX_EE_F="${MAX_EE_F:-2}"                # Max Expected Errors allowed for Forward read
export MAX_EE_R="${MAX_EE_R:-5}"                # Max Expected Errors allowed for Reverse read (relaxed for PE300)

# -------------------------------------------------------------------------
# Print Orchestrator Banner
# -------------------------------------------------------------------------
echo "========================================================================="
echo "        __  ___               __            _   __ _____                 "
echo "       /  |/  /____ _ _____  / /_ ___   __ / | / // ___/                 "
echo "      / /|_/ // __ \`// ___/ / __// _ \\ /_//  |/ / \\__ \\                  "
echo "     / /  / // /_/ // /__  / /_ /  __/   / /|  / ___/ /                  "
echo "    /_/  /_/ \\__,_/ \\___/  \\__/ \\___/   /_/ |_/ /____/                   "
echo "                                                                         "
echo "         Modular 16S rRNA Longitudinal Pipeline Orchestrator             "
echo "========================================================================="
echo "CPU Threads:            ${THREADS}"
echo "Raw Read Directory:     ${DATA_DIR}"
echo "Sample Sheet TSV:       ${SAMPLE_TSV}"
echo "Standard Metadata TSV:  ${METADATA_PATH}"
echo "V3-V4 Primer Forward:   ${PRIMER_F}"
echo "V3-V4 Primer Reverse:   ${PRIMER_R}"
echo "DADA2 Truncation F/R:   ${TRUNC_LEN_F} / ${TRUNC_LEN_R} bp"
echo "DADA2 Max EE F/R:       ${MAX_EE_F} / ${MAX_EE_R}"
echo "========================================================================="

# Helper Function: Display stage execution header
run_stage() {
    local stage_num="$1"
    local stage_name="$2"
    local stage_cmd="$3"
    
    echo ""
    echo "========================================================================="
    echo " >>> STAGE ${stage_num}: ${stage_name}"
    echo "     Running command: ${stage_cmd}"
    echo "========================================================================="
    
    # Execute command and track execution time
    local start_time=$(date +%s)
    
    if eval "${stage_cmd}"; then
        local end_time=$(date +%s)
        local elapsed=$((end_time - start_time))
        echo " >>> STAGE ${stage_num} completed successfully in ${elapsed}s."
    else
        echo "❌ ERROR: STAGE ${stage_num} (${stage_name}) failed!"
        exit 1
    fi
}

# Validate or auto-generate the sample TSV sheet
validate_or_generate_samples() {
    if [ ! -f "${SAMPLE_TSV}" ]; then
        echo ""
        echo "⚠️  WARNING: Standardized sample sheet not found at: ${SAMPLE_TSV}"
        echo ">>> Auto-generating sample sheet from FASTQ files in data folder..."
        bash script/generate_samples_tsv.sh "${DATA_DIR}" "${SAMPLE_TSV}"
    else
        echo "✅ Verified sample sheet exists at: ${SAMPLE_TSV}"
    fi
}

# Validate or auto-generate the standardized metadata file
validate_or_generate_metadata() {
    if [ ! -f "${METADATA_PATH}" ]; then
        echo ""
        echo "⚠️  WARNING: Standardized metadata file not found at: ${METADATA_PATH}"
        echo ">>> Automatically generating longitudinal metadata TSV based on sample sheet..."
        Rscript script/convert_metadata.R
    else
        echo "✅ Verified longitudinal metadata exists at: ${METADATA_PATH}"
    fi
}

# -------------------------------------------------------------------------
# 2. Pipeline Execution Sequence
# -------------------------------------------------------------------------

# Stage 0: Validate or generate input sample list
validate_or_generate_samples

# Stage 0.5: MD5 Raw Data Integrity Check
run_stage "0.5" "MD5 Raw Data Integrity Check" "bash script/check_md5.sh ${DATA_DIR}"

# Stage 1: Quality Control & Primer Trimming (Cutadapt, FastQC, MultiQC)
run_stage "1" "QC & Cutadapt Primer Trimming" "bash script/01_qc_trim.sh"

# Stage 2: Download Reference Databases (SILVA v138.1)
run_stage "2" "Reference Database Downloader" "bash script/download_db.sh"

# Stage 3: ASV Denoising & Taxonomy Assignment (DADA2 + SILVA)
run_stage "3" "DADA2 ASV Inference & Classification" "Rscript script/02_dada2.R"

# Stage 3.5: Validate or auto-generate metadata mapping based on processed amplicons
validate_or_generate_metadata

# Stage 4: Phyloseq Building & Table Export
run_stage "4" "Phyloseq Assembler & Excel Abundance Export" "Rscript script/03_phyloseq_prep.R"

# Stage 5: Longitudinal Repeated-Measures Statistics (LMM, PCoA, Volatility)
run_stage "5" "Repeated-Measures LMM & Beta Ordinations" "Rscript script/04_longitudinal_stats.R"

# Stage 6: Group Ecological Co-occurrence Networks
run_stage "6" "Ecological Interaction Network Mapping" "Rscript script/05_network_analysis.R"

# Stage 7: Metagenome Functional Prediction (PICRUSt2)
run_stage "7" "PICRUSt2 Metagenomic Functional Forecasting" "bash script/06_picrust2.sh"

echo ""
echo "========================================================================="
echo "🎉 SUCCESS: The entire Longitudinal 16S Workflow has finished!"
echo "All statistics, ecological networks, abundance sheets, and functional"
echo "predictions are compiled inside: ${RESULTS_DIR}/"
echo "========================================================================="
