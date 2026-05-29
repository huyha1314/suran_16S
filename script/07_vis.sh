#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define path variables
BASE_DIR="/worker_data2/huyha/precisiongene/suran_16S"
VIS_DIR="$BASE_DIR/vis"
RESULTS_DIR="$BASE_DIR/results"
RP_DIR="$BASE_DIR/rp_final"
INTERACTIVE_OUT="$RESULTS_DIR/interactive_reports"

echo "=========================================================="
echo " Starting 16S Report Generation & Environment Setup       "
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
cp "$RESULTS_DIR/04_stats/alpha/kruskal_shannon.tsv" "$RP_DIR/04_stats/alpha/"
cp "$RESULTS_DIR/04_stats/alpha/pairwise_wilcoxon_shannon.tsv" "$RP_DIR/04_stats/alpha/"

# Copy all Beta Diversity files (permanova, permdisp, and ordination PNGs)
cp "$RESULTS_DIR"/04_stats/beta/* "$RP_DIR/04_stats/beta/"

cp "$RESULTS_DIR"/04_stats/differential/* "$RP_DIR/04_stats/differential/"
cp "$RESULTS_DIR/06_picrust2/pathways_out/path_abun_unstrat.tsv" "$RP_DIR/06_picrust2/pathways_out/"

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

cp ./7.rp.qmd ./rp_final
cd rp_final
pixi run quarto render ./7.rp.qmd 
rm ./7.rp.qmd 
tar -czf ../16S_suran.tar.gz ./*