#!/usr/bin/env bash

# Exit immediately if any command exits with a non-zero status
set -euo pipefail

# Parameters
DATA_DIR="${1:-./data}"
OUTPUT_TSV="${2:-${DATA_DIR}/sample.tsv}"

echo "========================================================================="
echo "               16S Pipeline Sample TSV Generator                         "
echo "========================================================================="
echo "Input Data Directory: ${DATA_DIR}"
echo "Output TSV Path:      ${OUTPUT_TSV}"
echo "========================================================================="

# Ensure directory exists
mkdir -p "$(dirname "${OUTPUT_TSV}")"

# Create a temporary file
TEMP_FILE=$(mktemp)

# Write headers to TSV
echo -e "SampleID\tForwardPath\tReversePath" > "${TEMP_FILE}"

count=0

# Discover forward reads (R1/1)
# Supporting standard _R1_, _R1., _1.fastq.gz, _1.fq.gz patterns
for r1_path in $(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | sort); do
    filename=$(basename "${r1_path}")
    
    # Extract the SampleID (first field split by '_')
    # Matches the user's -d = "_" parameter
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
echo "========================================================================="
