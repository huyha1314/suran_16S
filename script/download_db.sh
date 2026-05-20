#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# Parameters
DB_DIR="${DB_DIR:-./data/db}"

echo "========================================================================="
echo "                16S Reference Database Downloader                        "
echo "========================================================================="
echo "Target Directory: ${DB_DIR}"
echo "========================================================================="

# Ensure target directory exists
mkdir -p "${DB_DIR}"

# SILVA Reference Database URLs (v138.1 formatted for DADA2)
SILVA_TRAIN_URL="https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1"
SILVA_SPECIES_URL="https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1"

# Target file paths
SILVA_TRAIN_FILE="${DB_DIR}/silva_nr99_v138.1_train_set.fa.gz"
SILVA_SPECIES_FILE="${DB_DIR}/silva_species_assignment_v138.1.fa.gz"

# Helper function to download files with resume support
download_db_file() {
    local url="$1"
    local dest="$2"
    local name="$3"
    
    if [ -f "${dest}" ] && [ -s "${dest}" ]; then
        echo "✅ ${name} already exists and is non-empty at: ${dest}"
    else
        echo "⏳ Downloading ${name} from Zenodo..."
        if command -v wget &> /dev/null; then
            # -c allows resuming broken downloads, -O sets output file
            wget -c -O "${dest}" "${url}"
        elif command -v curl &> /dev/null; then
            # -L follows redirects, -C - resumes download, -o sets output
            curl -L -C - -o "${dest}" "${url}"
        else
            echo "❌ ERROR: Neither 'wget' nor 'curl' was found in your system path."
            echo "Please install one of them to download database reference files."
            exit 1
        fi
        echo "✅ Successfully downloaded ${name}!"
    fi
}

# Execute downloads in parallel using background processes
download_db_file "${SILVA_TRAIN_URL}" "${SILVA_TRAIN_FILE}" "SILVA v138.1 Training Set" &
pid1=$!
download_db_file "${SILVA_SPECIES_URL}" "${SILVA_SPECIES_FILE}" "SILVA v138.1 Species Assignment" &
pid2=$!

# Wait for both processes to finish
wait $pid1
wait $pid2

echo "========================================================================="
echo "🎉 Reference database download and verification complete!"
echo "All databases are ready in: ${DB_DIR}"
echo "========================================================================="
