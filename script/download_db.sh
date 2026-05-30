#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# Parameters
DB_DIR="${DB_DIR:-./data/db}"
MODE="${MODE:-shortread}"

echo "========================================================================="
echo "                16S Reference Database Downloader                        "
echo "========================================================================="
echo "Target Directory: ${DB_DIR}"
echo "Execution Mode:   ${MODE}"
echo "========================================================================="

# Ensure target directory exists
mkdir -p "${DB_DIR}"

# SILVA Reference Database URLs (v138.1 formatted for DADA2)
SILVA_TRAIN_URL="https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1"
SILVA_SPECIES_URL="https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1"

# Target file paths
SILVA_TRAIN_FILE="${SILVA_TRAIN_FILE:-${DB_DIR}/silva_nr99_v138.1_train_set.fa.gz}"
SILVA_SPECIES_FILE="${SILVA_SPECIES_FILE:-${DB_DIR}/silva_species_assignment_v138.1.fa.gz}"

# Helper function to download files with resume support
download_db_file() {
    local url="$1"
    local dest="$2"
    local name="$3"
    
    if [ -f "${dest}" ] && [ -s "${dest}" ]; then
        echo "✅ ${name} already exists and is non-empty at: ${dest}"
    else
        # Ensure the destination directory exists
        mkdir -p "$(dirname "${dest}")"
        echo "⏳ Downloading ${name}..."

        if command -v wget &> /dev/null; then
            wget -c -O "${dest}" "${url}"
        elif command -v curl &> /dev/null; then
            curl -L -C - -o "${dest}" "${url}"
        else
            echo "❌ ERROR: Neither 'wget' nor 'curl' was found in your system path."
            echo "Please install one of them to download database reference files."
            exit 1
        fi
        echo "✅ Successfully downloaded ${name}!"
    fi
}

# 1. Download short-read/DADA2 SILVA databases
download_db_file "${SILVA_TRAIN_URL}" "${SILVA_TRAIN_FILE}" "SILVA v138.1 Training Set" &
pid1=$!
download_db_file "${SILVA_SPECIES_URL}" "${SILVA_SPECIES_FILE}" "SILVA v138.1 Species Assignment" &
pid2=$!

# Wait for SILVA downloads
wait $pid1
wait $pid2

# 2. If longread mode, also download Emu and Kraken2/Bracken databases
if [ "${MODE}" == "longread" ]; then
    echo ""
    echo "------------------- Long-read Databases -------------------"
    
    # Emu Database Setup
    EMU_DIR="${DB_DIR}/emu"
    if [ -f "${EMU_DIR}/species_taxid.fasta" ] && [ -f "${EMU_DIR}/taxonomy.tsv" ]; then
        echo "✅ Emu default database already exists at: ${EMU_DIR}"
    else
        echo "⏳ Downloading Emu default database using osfclient..."
        mkdir -p "${EMU_DIR}"
        (
            cd "${EMU_DIR}"
            pixi run osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
            echo "📦 Extracting Emu database..."
            tar -xf emu.tar
            rm -f emu.tar
        )
        echo "✅ Emu database setup completed successfully!"
    fi

    # Kraken2/Bracken SILVA 16S Database Setup
    KRAKEN_DIR="${DB_DIR}/kraken2"
    if [ -f "${KRAKEN_DIR}/hash.k2d" ] && [ -f "${KRAKEN_DIR}/opts.k2d" ]; then
        echo "✅ Kraken2/Bracken SILVA 16S database already exists at: ${KRAKEN_DIR}"
    else
        echo "⏳ Downloading Kraken2/Bracken SILVA 16S database..."
        mkdir -p "${KRAKEN_DIR}"
        (
            cd "${KRAKEN_DIR}"
            if command -v wget &> /dev/null; then
                wget -c https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz
            else
                curl -L -O https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz
            fi
            echo "📦 Extracting Kraken2/Bracken database..."
            tar -xf 16S_Silva138_20200326.tgz
            rm -f 16S_Silva138_20200326.tgz
        )
        echo "✅ Kraken2/Bracken database setup completed successfully!"
    fi
fi

echo "========================================================================="
echo "🎉 Reference database download and verification complete!"
echo "All databases are ready in: ${DB_DIR}"
echo "========================================================================="
