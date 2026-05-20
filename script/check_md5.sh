#!/usr/bin/env bash

# Exit immediately if any command exits with a non-zero status
set -euo pipefail

# Parameters
DATA_DIR="${1:-./data}"

echo "========================================================================="
echo "               16S Raw Sequencing MD5 Verification Tool                 "
echo "========================================================================="
echo "Scanning Directory: ${DATA_DIR}"
echo "========================================================================="

if [ ! -d "${DATA_DIR}" ]; then
    echo "ERROR: Data directory not found: ${DATA_DIR}"
    exit 1
fi

# Find MD5 signature files
md5_files=$(find "${DATA_DIR}" -maxdepth 1 -name "*.md5" | sort)

if [ -z "${md5_files}" ]; then
    echo "⚠️  WARNING: No .md5 hash files found in ${DATA_DIR}."
    echo "Running fallback integrity check (testing gzip compression)..."
    
    fastq_files=$(find "${DATA_DIR}" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
    if [ -z "${fastq_files}" ]; then
        echo "❌ ERROR: No compressed FASTQ files (*.fastq.gz or *.fq.gz) found in ${DATA_DIR}."
        exit 1
    fi
    
    failed_gzip=0
    for fq in ${fastq_files}; do
        echo -n "Checking Gzip Integrity: $(basename "${fq}") ... "
        if gzip -t "${fq}" 2>/dev/null; then
            echo "OK"
        else
            echo "FAILED (corrupted or truncated file!)"
            failed_gzip=$((failed_gzip + 1))
        fi
    done
    
    echo "-------------------------------------------------------------------------"
    if [ "${failed_gzip}" -gt 0 ]; then
        echo "❌ MD5 FALLBACK FAILURE: ${failed_gzip} file(s) failed compression check!"
        exit 1
    else
        echo "🎉 MD5 FALLBACK SUCCESS: All compressed reads are valid gzip archives!"
        exit 0
    fi
fi

# MD5 list exists
file_count=$(echo "${md5_files}" | wc -w)
echo "Found ${file_count} MD5 signature files. Starting verification..."
echo "-------------------------------------------------------------------------"

failed_count=0
success_count=0

# Save current directory to return to it later
orig_dir=$(pwd)

# Run inside the target directory so relative paths in md5 files resolve correctly
cd "${DATA_DIR}"

for md5_file in $(find . -maxdepth 1 -name "*.md5" | sort); do
    filename=$(basename "${md5_file}")
    
    # Read the signature line and strip Windows/DOS carriage returns (\r)
    signature_line=$(tr -d '\r' < "${filename}" | xargs)
    
    if [ -z "${signature_line}" ]; then
        echo "⚠️  $(basename "${md5_file}"): MD5 file is empty, skipping..."
        continue
    fi
    
    # Extract expected hash (first field)
    expected_hash=$(echo "${signature_line}" | awk '{print $1}')
    
    # Extract expected target filename (everything after the first field, cleaned)
    target_filename=$(echo "${signature_line}" | awk '{$1=""; print $0}' | sed -e 's/^[ \t]*//')
    
    if [ -z "${target_filename}" ]; then
        # Fallback if parsing space is non-standard
        target_filename="${filename%.md5}"
    fi
    
    if [ -f "${target_filename}" ]; then
        # Calculate actual md5 sum
        actual_hash=$(md5sum "${target_filename}" | awk '{print $1}')
        
        if [ "${expected_hash}" = "${actual_hash}" ]; then
            echo "✅ ${target_filename} : OK"
            success_count=$((success_count + 1))
        else
            echo "❌ ${target_filename} : MD5 MISMATCH!"
            echo "   Expected: ${expected_hash}"
            echo "   Actual:   ${actual_hash}"
            failed_count=$((failed_count + 1))
        fi
    else
        echo "❌ ${target_filename} : FILE NOT FOUND!"
        failed_count=$((failed_count + 1))
    fi
done

# Return to starting directory
cd "${orig_dir}"

echo "-------------------------------------------------------------------------"
echo "Verification Summary:"
echo "  Total Verified:  $((success_count + failed_count))"
echo "  Success:         ${success_count}"
echo "  Failed:          ${failed_count}"
echo "========================================================================="

if [ "${failed_count}" -gt 0 ]; then
    echo "❌ ERROR: One or more sequencing files are corrupted or incomplete!"
    exit 1
else
    echo "🎉 SUCCESS: All raw FASTQ files verified successfully against MD5 hashes!"
fi
