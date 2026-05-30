#!/bin/bash
# rp.sh — Quick report rebuild & packaging script
# All settings are read from pipeline.config (same source as master.sh).
# Usage:  bash rp.sh
# Or override a single variable: MODE=shortread bash rp.sh

set -e

BASE_DIR="$PWD"
RP_DIR="$BASE_DIR/rp_final"

# ── Load all settings from pipeline.config ────────────────────────────────
if [ -f "$BASE_DIR/pipeline.config" ]; then
    source "$BASE_DIR/pipeline.config"
else
    echo "Warning: pipeline.config not found — using built-in defaults."
fi

# Fallback defaults (mirror master.sh)
export MODE="${MODE:-shortread}"
export PROJECT_NAME="${PROJECT_NAME:-16S_analysis}"
export RESULTS_DIR="${RESULTS_DIR:-./results}"

# ── Select QMD template based on MODE ─────────────────────────────────────
if [ "${MODE}" == "longread" ]; then
    QMD_TEMPLATE="$BASE_DIR/7.rp_longread.qmd"
    echo "-> Report mode: Long-Read ONT (Emu + Bracken)"
else
    QMD_TEMPLATE="$BASE_DIR/7.rp_shortread.qmd"
    echo "-> Report mode: Short-Read Illumina (DADA2 + PICRUSt2)"
fi

if [ ! -f "${QMD_TEMPLATE}" ]; then
    echo "ERROR: QMD template not found at ${QMD_TEMPLATE}"
    exit 1
fi

# ── Copy template and render ───────────────────────────────────────────────
echo "-> Copying template to rp_final/..."
cp "${QMD_TEMPLATE}" "${RP_DIR}/7.rp.qmd"

echo "-> Rendering Quarto report..."
cd "${RP_DIR}"
pixi run quarto render ./7.rp.qmd
rm ./7.rp.qmd

# ── Package into archives ──────────────────────────────────────────────────
echo "-> Packaging report as archive..."
cd "${BASE_DIR}"
tar -czf "${PROJECT_NAME}.tar.gz" -C rp_final .
zip -r "${PROJECT_NAME}.zip" rp_final/ -x "*.DS_Store" -x "*/__pycache__/*"

echo "=========================================================="
echo " Report built successfully!"
echo " HTML report:  rp_final/7.rp.html"
echo " Archive (tar): ${PROJECT_NAME}.tar.gz"
echo " Archive (zip): ${PROJECT_NAME}.zip"
echo "=========================================================="