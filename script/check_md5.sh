#!/usr/bin/env bash

# Exit immediately if any command exits with a non-zero status
set -euo pipefail

# Parameters
DATA_DIR="${1:-./data}"
SAMPLE_TSV="${SAMPLE_TSV:-./sample.tsv}"
THREADS="${THREADS:-8}"
MARKER_FILE="log/.md5_verified"

echo "========================================================================="
# Check if GNU parallel is installed
if command -v parallel &> /dev/null; then
    PARALLEL_TOOL="GNU Parallel"
else
    PARALLEL_TOOL="xargs (fallback)"
fi

echo "        16S Raw Sequencing Verification Tool (Strict Sample TSV)        "
echo "========================================================================="
echo "Sample Sheet TSV:  ${SAMPLE_TSV}"
echo "Threads Allocated: ${THREADS}"
echo "Parallel Tool:     ${PARALLEL_TOOL}"
echo "========================================================================="

# Validate sample.tsv exists
if [ ! -f "${SAMPLE_TSV}" ]; then
    echo "❌ ERROR: Standardized sample sheet not found at: ${SAMPLE_TSV}"
    echo "Please generate it first by running: pixi run generate_samples"
    exit 1
fi

# Extract unique raw fastq paths from columns 2 and 3 of the sample sheet
files_to_check=$(tail -n +2 "${SAMPLE_TSV}" | cut -f2,3 | tr '\t' '\n' | sort -u | grep -v '^$')

if [ -z "${files_to_check}" ]; then
    echo "❌ ERROR: No raw FASTQ file paths found in ${SAMPLE_TSV}."
    exit 1
fi

file_count=$(echo "${files_to_check}" | wc -w)
echo "Found ${file_count} unique files listed in sample sheet. Starting parallel verification..."
echo "-------------------------------------------------------------------------"

# Define the parallel verification worker function
verify_file() {
    local fq="$1"
    
    if [ ! -f "${fq}" ]; then
        echo "❌ ${fq} : FILE NOT FOUND!"
        return 1
    fi
    
    local md5_file="${fq}.md5"
    if [ -f "${md5_file}" ]; then
        # Run parallel-safe MD5 check (strip carriage returns first)
        if tr -d '\r' < "${md5_file}" | (cd "$(dirname "${fq}")" && md5sum -c - &>/dev/null); then
            echo "✅ $(basename "${fq}") : MD5 OK"
            return 0
        else
            echo "❌ $(basename "${fq}") : MD5 MISMATCH!"
            return 1
        fi
    else
        # Run fallback Gzip integrity check
        if gzip -t "${fq}" 2>/dev/null; then
            echo "✅ $(basename "${fq}") : GZIP OK (No MD5 file)"
            return 0
        else
            echo "❌ $(basename "${fq}") : GZIP CORRUPTED!"
            return 1
        fi
    fi
}
export -f verify_file

# Run parallel verification
verification_failed=0

if command -v parallel &> /dev/null; then
    # GNU Parallel approach
    if parallel -j "${THREADS}" verify_file ::: ${files_to_check}; then
        echo "-------------------------------------------------------------------------"
        echo "🎉 SUCCESS: All raw FASTQ files verified successfully!"
    else
        echo "-------------------------------------------------------------------------"
        echo "❌ ERROR: One or more sequencing files are missing, corrupted, or have MD5 mismatches!"
        verification_failed=1
    fi
else
    # xargs fallback approach
    if echo "${files_to_check}" | xargs -I {} -P "${THREADS}" bash -c 'verify_file "{}"'; then
        echo "-------------------------------------------------------------------------"
        echo "🎉 SUCCESS: All raw FASTQ files verified successfully!"
    else
        echo "-------------------------------------------------------------------------"
        echo "❌ ERROR: One or more sequencing files are missing, corrupted, or have MD5 mismatches!"
        verification_failed=1
    fi
fi

if [ "${verification_failed}" -eq 0 ]; then
    # Create the verification marker so it is cached for future runs
    mkdir -p "$(dirname "${MARKER_FILE}")"
    echo "verified $(date)" > "${MARKER_FILE}"
    echo "🔒 Verification marker created at: ${MARKER_FILE}"
    exit 0
else
    exit 1
fi
