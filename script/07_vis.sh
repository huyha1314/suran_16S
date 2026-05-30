#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define path variables
BASE_DIR="$PWD"

# ── Load all settings from pipeline.config ────────────────────────────────
if [ -f "$BASE_DIR/pipeline.config" ]; then
    source "$BASE_DIR/pipeline.config"
fi

# Fallback defaults
export MODE="${MODE:-shortread}"
export RESULTS_DIR="${RESULTS_DIR:-$BASE_DIR/results}"
export PROJECT_NAME="${PROJECT_NAME:-16S_analysis}"

VIS_DIR="$BASE_DIR/vis"
RP_DIR="$BASE_DIR/rp_final"
INTERACTIVE_OUT="$RESULTS_DIR/interactive_reports"

echo "=========================================================="
echo " Starting 16S Report Generation & Environment Setup       "
echo " Mode: ${MODE} | Project: ${PROJECT_NAME}"
echo "=========================================================="
[ -d "$INTERACTIVE_OUT" ] && rm -r "$INTERACTIVE_OUT"
# Step 1: Ensure directory structure exists
echo "-> Creating necessary output directories..."
mkdir -p "$INTERACTIVE_OUT"
mkdir -p "$RP_DIR/interactive_reports"
mkdir -p "$RP_DIR/03_phyloseq"
mkdir -p "$RP_DIR/04_stats/alpha"
mkdir -p "$RP_DIR/04_stats/beta"
mkdir -p "$RP_DIR/04_stats/differential"
mkdir -p "$RP_DIR/06_picrust2/pathways_out"

# Step 2: Run all visualization R scripts sequentially
echo "-> Running standalone R visualization scripts..."

echo "   Running 1.tax.R..."
Rscript "$VIS_DIR/1.tax.R"

echo "   Running 2.alpha.R..."
Rscript "$VIS_DIR/2.alpha.R"

echo "   Running 3.beta.R..."
Rscript "$VIS_DIR/3.beta.R"

echo "   Running 4.rarefaction.R..."
Rscript "$VIS_DIR/4.rarefaction.R"

echo "   Running 5.differential.R..."
Rscript "$VIS_DIR/5.differential.R"

echo "   Running 6.function.R..."
Rscript "$VIS_DIR/6.function.R"

echo "-> All R visualizations generated successfully."

# Step 3: Copy generated HTML widgets and high-res PDF plots to the final report folder
echo "-> Copying interactive HTML files and PDF plots to the report directory..."
cp -r $INTERACTIVE_OUT/* "$RP_DIR/interactive_reports/"

# Step 4: Copy raw data inputs for the download buttons
echo "-> Syncing download data files to the report directory..."
cp "$RESULTS_DIR/multiqc_report.html" "$RP_DIR/"
cp "$RESULTS_DIR/03_phyloseq/taxa_abundance_tables.xlsx" "$RP_DIR/03_phyloseq/"
cp "$RESULTS_DIR"/03_phyloseq/abundance_*.tsv "$RP_DIR/03_phyloseq/"
cp "$RESULTS_DIR/04_stats/alpha/alpha_diversity_table.tsv" "$RP_DIR/04_stats/alpha/"
[ -f "$RESULTS_DIR/04_stats/alpha/kruskal_shannon.tsv" ] && cp "$RESULTS_DIR/04_stats/alpha/kruskal_shannon.tsv" "$RP_DIR/04_stats/alpha/"
[ -f "$RESULTS_DIR/04_stats/alpha/pairwise_wilcoxon_shannon.tsv" ] && cp "$RESULTS_DIR/04_stats/alpha/pairwise_wilcoxon_shannon.tsv" "$RP_DIR/04_stats/alpha/"

# Copy all Beta Diversity files (permanova, permdisp, and ordination PNGs)
cp "$RESULTS_DIR"/04_stats/beta/* "$RP_DIR/04_stats/beta/"

cp "$RESULTS_DIR"/04_stats/differential/* "$RP_DIR/04_stats/differential/"

# Copy PICRUSt2 output only if it exists (skipped in longread mode)
if [ -f "$RESULTS_DIR/06_picrust2/pathways_out/path_abun_unstrat.tsv" ]; then
    mkdir -p "$RP_DIR/06_picrust2/pathways_out/"
    cp "$RESULTS_DIR/06_picrust2/pathways_out/path_abun_unstrat.tsv" "$RP_DIR/06_picrust2/pathways_out/"
else
    echo "   Note: PICRUSt2 output not found (expected in longread mode). Skipping."
fi

# Step 5: Copy UI headers/footers
echo "-> Copying header and footer styling folders..."
if [ -d "$VIS_DIR/header_footer" ]; then
    cp -r "$VIS_DIR/header_footer" "$RP_DIR/"
else
    echo "Warning: header_footer directory not found in $VIS_DIR"
fi

echo "=========================================================="
echo " Environment setup complete!                              "
echo " Move your 16S_report.qmd into: $RP_DIR/                  "
echo " Run: 'quarto render 16S_report.qmd' inside that folder.  "
echo "=========================================================="

# Step 6: Select the correct Quarto report template based on pipeline mode
if [ "${MODE}" == "longread" ]; then
    QMD_TEMPLATE="./7.rp_longread.qmd"
    REPORT_NAME="7.rp.qmd"
    echo "-> Selected report template: Long-Read ONT (${QMD_TEMPLATE})"
else
    QMD_TEMPLATE="./7.rp_shortread.qmd"
    REPORT_NAME="7.rp.qmd"
    echo "-> Selected report template: Short-Read Illumina (${QMD_TEMPLATE})"
fi

cp "${QMD_TEMPLATE}" "${RP_DIR}/${REPORT_NAME}"
cd rp_final
pixi run quarto render "./${REPORT_NAME}"
rm "./${REPORT_NAME}"
tar -czf "../${PROJECT_NAME}.tar.gz" ./*
zip -r "../${PROJECT_NAME}.zip" ./* -x "*.DS_Store"

echo "=========================================================="
echo " Report build complete!"
echo " HTML:    rp_final/7.rp.html"
echo " Archive: ${PROJECT_NAME}.tar.gz / ${PROJECT_NAME}.zip"
echo "=========================================================="