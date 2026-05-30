#!/usr/bin/env bash

# Exit immediately if any command exits with a non-zero status
set -euo pipefail

# Parameters
DATA_DIR="${1:-./data}"
OUTPUT_TSV="${2:-${DATA_DIR}/sample.tsv}"

# Load external configuration if present to discover pipeline MODE, while preserving environment overrides
PRESERVED_MODE="${MODE:-}"
if [ -f "pipeline.config" ]; then
    source pipeline.config
fi
if [ -n "${PRESERVED_MODE}" ]; then
    export MODE="${PRESERVED_MODE}"
fi
export MODE="${MODE:-shortread}"

echo "========================================================================="
echo "               16S Pipeline Sample TSV Generator                         "
echo "========================================================================="
echo "Input Data Directory: ${DATA_DIR}"
echo "Output TSV Path:      ${OUTPUT_TSV}"
echo "Pipeline Mode:        ${MODE}"
echo "========================================================================="

# Ensure directory exists
mkdir -p "$(dirname "${OUTPUT_TSV}")"

# Create a temporary file
TEMP_FILE=$(mktemp)

# Write headers to TSV
echo -e "SampleID\tForwardPath\tReversePath" > "${TEMP_FILE}"

count=0

if [ "${MODE}" == "longread" ]; then
    # Discover single-end long reads (ONT / PacBio)
    # Supporting any .fastq.gz or .fq.gz files
    for r1_path in $(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort); do
        filename=$(basename "${r1_path}")
        
        # Skip sample.tsv just in case
        [ "${filename}" == "sample.tsv" ] && continue
        
        # Extract the SampleID by stripping common long-read suffixes
        # This preserves compound names with underscores (e.g., DC3_01, EM_10)
        sample_id="${filename}"
        for suffix in "_trimmed_1.fastq.gz" "_trimmed_1.fq.gz" "_trimmed.fastq.gz" "_trimmed.fq.gz" "_1.fastq.gz" "_1.fq.gz" ".fastq.gz" ".fq.gz"; do
            if [[ "${sample_id}" == *"${suffix}" ]]; then
                sample_id="${sample_id%${suffix}}"
                break
            fi
        done
        
        # Write single-end forward read with blank reverse read path
        echo -e "${sample_id}\t${r1_path}\t" >> "${TEMP_FILE}"
        count=$((count + 1))
    done
    
    # Save final TSV
    mv "${TEMP_FILE}" "${OUTPUT_TSV}"
    echo "SUCCESS: Generated sample TSV at ${OUTPUT_TSV} with ${count} single-end long-read samples!"
else
    # Discover paired-end short reads (Illumina V3-V4)
    for r1_path in $(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | sort); do
        filename=$(basename "${r1_path}")
        
        # Extract the SampleID (first field split by '_')
        sample_id=$(echo "${filename}" | cut -d'_' -f1)
        
        # Construct corresponding paired-end R2 path
        r2_path=""
        if [[ "${filename}" == *"_R1_"* ]]; then
            r2_path="${r1_path/_R1_/_R2_}"
        elif [[ "${filename}" == *"_R1."* ]]; then
            r2_path="${r1_path/_R1./_R2.}"
        elif [[ "${filename}" == *"_1.fastq.gz"* ]]; then
            r2_path="${r1_path/_1.fastq.gz/_2.fastq.gz}"
        elif [[ "${filename}" == *"_1.fq.gz"* ]]; then
            r2_path="${r1_path/_1.fq.gz/_2.fq.gz}"
        fi
        
        # Check if the paired reverse read file actually exists
        if [ -f "${r2_path}" ]; then
            echo -e "${sample_id}\t${r1_path}\t${r2_path}" >> "${TEMP_FILE}"
            count=$((count + 1))
        else
            echo "Warning: No matching reverse read found at: ${r2_path} for: ${filename}"
        fi
    done
    
    # Save final TSV
    mv "${TEMP_FILE}" "${OUTPUT_TSV}"
    echo "SUCCESS: Generated sample TSV at ${OUTPUT_TSV} with ${count} paired samples!"
fi

echo "========================================================================="

