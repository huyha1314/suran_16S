#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status, 
# if an uninitialized variable is used, or if a piped command fails.
set -euo pipefail

# Setup automated run-level logging inside log/ directory
LOG_DIR="log"
mkdir -p "${LOG_DIR}"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/pipeline_${TIMESTAMP}.log"

# Redirect stdout and stderr to both console and the timestamped log file
exec > >(tee -i "${LOG_FILE}") 2>&1

echo "========================================================================="
echo "📝 LOGGING ACTIVE: Live output is being saved to: ${LOG_FILE}"
echo "========================================================================="

# =========================================================================
#      LONGITUDINAL 16S WORKFLOW ORCHESTRATOR & CONFIGURATION PANEL
# =========================================================================

# Load external configuration if present
if [ -f "pipeline.config" ]; then
    echo "⚙️ Loading external configuration from pipeline.config..."
    source pipeline.config
fi

# -------------------------------------------------------------------------
# Fallback defaults (pipeline.config values take priority when sourced above)
# -------------------------------------------------------------------------
# Project identity
export PROJECT_NAME="${PROJECT_NAME:-16S_analysis}"

# Threads & concurrency
export THREADS="${THREADS:-64}"
export MD5_THREADS="${MD5_THREADS:-${THREADS}}"
export FASTQC_THREADS="${FASTQC_THREADS:-${THREADS}}"
export TRIM_CONCURRENCY="${TRIM_CONCURRENCY:-8}"
export TRIM_THREADS_PER_JOB="${TRIM_THREADS_PER_JOB:-8}"
export DADA2_THREADS="${DADA2_THREADS:-${THREADS}}"
export PICRUST2_THREADS="${PICRUST2_THREADS:-${THREADS}}"

# Input paths
export DATA_DIR="${DATA_DIR:-./data}"
export SAMPLE_TSV="${SAMPLE_TSV:-./sample.tsv}"
export METADATA_PATH="${METADATA_PATH:-./metadata.tsv}"
export DB_DIR="${DB_DIR:-./data/db}"
export SILVA_TRAIN_FILE="${SILVA_TRAIN_FILE:-${DB_DIR}/silva_nr99_v138.1_train_set.fa.gz}"
export SILVA_SPECIES_FILE="${SILVA_SPECIES_FILE:-${DB_DIR}/silva_species_assignment_v138.1.fa.gz}"

# Output directories
export RESULTS_DIR="${RESULTS_DIR:-./results}"
export QC_DIR="${QC_DIR:-${RESULTS_DIR}/01_trimmed/qc}"
export TRIM_DIR="${TRIM_DIR:-${RESULTS_DIR}/01_trimmed}"
export DADA2_DIR="${DADA2_DIR:-${RESULTS_DIR}/02_dada2}"
export PHYLO_DIR="${PHYLO_DIR:-${RESULTS_DIR}/03_phyloseq}"
export STATS_DIR="${STATS_DIR:-${RESULTS_DIR}/04_stats}"
export NETWORK_DIR="${NETWORK_DIR:-${RESULTS_DIR}/05_network}"
export PICRUSt2_DIR="${PICRUSt2_DIR:-${RESULTS_DIR}/06_picrust2}"

# Pipeline mode
export MODE="${MODE:-shortread}"

# Primer sequences
export PRIMER_F="${PRIMER_F:-CCTACGGGNGGCWGCAG}"
export PRIMER_R="${PRIMER_R:-GACTACHVGGGTATCTAATCC}"
export PRIMER_F_LONG="${PRIMER_F_LONG:-AGAGTTTGATCMTGGCTCAG}"
export PRIMER_R_LONG="${PRIMER_R_LONG:-CGGTTACCTTGTTACGACTT}"

# DADA2 denoising parameters
export TRUNC_LEN_F="${TRUNC_LEN_F:-270}"
export TRUNC_LEN_R="${TRUNC_LEN_R:-210}"
export MAX_EE_F="${MAX_EE_F:-2}"
export MAX_EE_R="${MAX_EE_R:-5}"

# Filtlong long-read quality filters
export FILTLONG_MIN_LEN="${FILTLONG_MIN_LEN:-800}"
export FILTLONG_MIN_Q="${FILTLONG_MIN_Q:-10}"


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
echo "Global Threads fallback: ${THREADS}"
echo "------------------- Parallel Resource Allocations ---------------------"
echo "1. MD5 Integrity Check: ${MD5_THREADS} cores (Parallel)"
echo "2. Raw/Trimmed FastQC:  ${FASTQC_THREADS} cores (Parallel)"
echo "3. Cutadapt Trimming:   ${TRIM_CONCURRENCY} concurrent jobs x ${TRIM_THREADS_PER_JOB} threads/job"
echo "4. DADA2 Denoising:     ${DADA2_THREADS} cores (Parallel)"
echo "5. PICRUSt2 Metagenome: ${PICRUST2_THREADS} cores (Parallel)"
echo "-------------------------------------------------------------------------"
echo "Raw Read Directory:     ${DATA_DIR}"
echo "Sample Sheet TSV:       ${SAMPLE_TSV}"
echo "Standard Metadata TSV:  ${METADATA_PATH}"
echo "Pipeline Execution Mode: ${MODE}"
if [ "${MODE}" == "longread" ]; then
    echo "ONT Primer Forward:     ${PRIMER_F_LONG}"
    echo "ONT Primer Reverse:     ${PRIMER_R_LONG}"
else
    echo "V3-V4 Primer Forward:   ${PRIMER_F}"
    echo "V3-V4 Primer Reverse:   ${PRIMER_R}"
    echo "DADA2 Truncation F/R:   ${TRUNC_LEN_F} / ${TRUNC_LEN_R} bp"
    echo "DADA2 Max EE F/R:       ${MAX_EE_F} / ${MAX_EE_R}"
fi
echo "========================================================================="

# Helper Function: Display stage execution header and support resume checks
run_stage() {
    local stage_num="$1"
    local stage_name="$2"
    local stage_cmd="$3"
    local check_file="${4:-}"
    
    if [ -n "${check_file}" ] && [ -f "${check_file}" ] && [ -s "${check_file}" ]; then
        echo ""
        echo "⏭️  [RESUME] STAGE ${stage_num}: ${stage_name} already completed."
        echo "    Verified output file exists: ${check_file}"
        echo "    Skipping to next stage..."
        return 0
    fi
    
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

# Validate that the standardized sample sheet exists
validate_samples_exist() {
    if [ ! -f "${SAMPLE_TSV}" ]; then
        echo ""
        echo "❌ ERROR: Standardized sample sheet not found at: ${SAMPLE_TSV}"
        echo "Please generate it manually first by running: pixi run generate_samples"
        exit 1
    else
        echo "✅ Verified sample sheet exists at: ${SAMPLE_TSV}"
    fi
}

# Validate that the standardized metadata file exists
validate_metadata_exist() {
    if [ ! -f "${METADATA_PATH}" ]; then
        echo ""
        echo "❌ ERROR: Standardized metadata file not found at: ${METADATA_PATH}"
        echo "Please generate it manually first by running: pixi run generate_metadata"
        exit 1
    else
        echo "✅ Verified longitudinal metadata exists at: ${METADATA_PATH}"
    fi
}

# -------------------------------------------------------------------------
# 2. Pipeline Execution Sequence
# -------------------------------------------------------------------------

# Stage 0: Strict Pre-flight Checks (Validate that manually generated TSVs exist)
validate_samples_exist
validate_metadata_exist

# Stage 0.5: MD5 Raw Data Integrity Check
run_stage "0.5" "MD5 Raw Data Integrity Check" "THREADS=${MD5_THREADS} bash script/check_md5.sh ${DATA_DIR}" "log/.md5_verified"

# Stage 1: Quality Control & Primer Trimming
if [ "${MODE}" == "longread" ]; then
    run_stage "1" "Long-read QC & Primer Trimming" "bash script/01_qc_trim_longread.sh" "${RESULTS_DIR}/multiqc_report.html"
else
    run_stage "1" "QC & Cutadapt Primer Trimming" "bash script/01_qc_trim.sh" "${RESULTS_DIR}/multiqc_report.html"
fi

# Stage 2: Download Reference Databases (SILVA v138.1)
# Note: Stage 2 remains identical as both modes map/train against the SILVA database
run_stage "2" "Reference Database Downloader" "bash script/download_db.sh" "${SILVA_SPECIES_FILE}"

# Stage 3: ASV Denoising / Long-read Denoising
if [ "${MODE}" == "longread" ]; then
    run_stage "3" "Long-read Taxonomy Profiling" "THREADS=${DADA2_THREADS} Rscript script/02_dada2_longread.R" "${DADA2_DIR}/seqtab_nochim.rds"
    run_stage "3.5" "Long-read Matrix Preparation (Emu & Kraken-Biom)" "mkdir -p ${RESULTS_DIR}/emu_combined && pixi run emu combine-outputs --counts ${DADA2_DIR}/emu species && cp ${DADA2_DIR}/emu/emu-combined-species-counts.tsv ${RESULTS_DIR}/emu_combined/emu-combined-abundance.tsv && pixi run kraken-biom ${DADA2_DIR}/kraken2/*_bracken.report --fmt json -o ${RESULTS_DIR}/bracken_combined.biom" "${RESULTS_DIR}/bracken_combined.biom"
else
    run_stage "3" "DADA2 ASV Inference & Classification" "THREADS=${DADA2_THREADS} Rscript script/02_dada2.R" "${DADA2_DIR}/seqtab_nochim.rds"
fi

# Stage 4: Phyloseq Building & Table Export
run_stage "4" "Phyloseq Assembler & Excel Abundance Export" "Rscript script/03_phyloseq_prep.R" "${PHYLO_DIR}/phyloseq_obj.rds"

# Stage 4.5: Decontam (Statistical Contaminant Removal)
run_stage "4.5" "Statistical Contaminant Removal (Decontam)" "Rscript script/03b_decontam_filter.R" ""

# Stage 5: Longitudinal Repeated-Measures Statistics (LMM, PCoA, Volatility)
run_stage "5" "Repeated-Measures LMM & Beta Ordinations" "Rscript script/04_longitudinal_stats.R" "${STATS_DIR}/differential/pairwise_differential_abundance.tsv"

# Stage 6: Group Ecological Co-occurrence Networks
run_stage "6" "Ecological Interaction Network Mapping" "Rscript script/05_network_analysis.R" "${NETWORK_DIR}/network_topology_comparison.tsv"

# Stage 7: Metagenome Functional Prediction (PICRUSt2)
if [ "${MODE}" != "longread" ]; then
    run_stage "7" "PICRUSt2 Metagenomic Functional Forecasting" "THREADS=${PICRUST2_THREADS} bash script/06_picrust2.sh" "${PICRUSt2_DIR}/pathways_out/path_abun_unstrat.tsv"
else
    echo "========================================================================="
    echo " >>> STAGE 7: Metagenome Functional Prediction (PICRUSt2)"
    echo "     SKIPPED: PICRUSt2 requires ASV sequences (not generated in longread mode)"
    echo "========================================================================="
fi

echo ""
echo "========================================================================="
echo "🎉 SUCCESS: The entire Longitudinal 16S Workflow has finished!"
echo "All statistics, ecological networks, abundance sheets, and functional"
echo "predictions are compiled inside: ${RESULTS_DIR}/"
echo "========================================================================="
