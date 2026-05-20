#!/usr/bin/env Rscript

# =========================================================================
#                 16S rRNA Phyloseq Object Generation & Export
# =========================================================================

cat("Loading required libraries...\n")
library(phyloseq)
library(tidyverse)
library(openxlsx)

# Get parameters from environment variables (with default fallbacks)
DADA2_DIR <- Sys.getenv("DADA2_DIR", "./results/02_dada2")
PHYLO_DIR <- Sys.getenv("PHYLO_DIR", "./results/03_phyloseq")
METADATA_PATH <- Sys.getenv("METADATA_PATH", "./data/metadata.tsv")

cat("=========================================================================\n")
cat("Parameters:\n")
cat("DADA2 Output Dir:   ", DADA2_DIR, "\n")
cat("Phyloseq Output Dir:", PHYLO_DIR, "\n")
cat("Metadata File Path: ", METADATA_PATH, "\n")
cat("=========================================================================\n")

# Create output directory
dir.create(PHYLO_DIR, showWarnings = FALSE, recursive = TRUE)

# Load DADA2 outputs
seqtab_path <- file.path(DADA2_DIR, "seqtab_nochim.rds")
taxa_path <- file.path(DADA2_DIR, "taxonomy.rds")

if (!file.exists(seqtab_path) || !file.exists(taxa_path)) {
  stop("ERROR: DADA2 output files not found. Please run DADA2 step first.")
}

seqtab <- readRDS(seqtab_path)
taxa <- readRDS(taxa_path)

sample_names_seqtab <- rownames(seqtab)

# -------------------------------------------------------------------------
# Step 1: Read Standardized Metadata TSV
# -------------------------------------------------------------------------
if (!file.exists(METADATA_PATH)) {
  cat("ERROR: Standardized metadata file not found at:", METADATA_PATH, "\n")
  cat("Please provide a tab-separated values (TSV) sheet with the following columns:\n")
  cat("  1. SampleID  - Name of the sample matching the raw FASTQ file name prefixes.\n")
  cat("  2. Group     - Treatment group (e.g., Group1, Group2, Group3).\n")
  cat("  3. Day       - Collection timepoint (e.g., Day1, Day2, Day3).\n")
  cat("  4. SubjectID - Unique ID of the subject to model repeated measures.\n")
  stop("FATAL: Standardized metadata file is missing!")
}

cat(">>> Reading standardized metadata TSV from:", METADATA_PATH, "\n")
metadata <- read_tsv(METADATA_PATH)

# Check required columns
required_cols <- c("SampleID", "Group", "Day", "SubjectID")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("ERROR: Standardized metadata is missing required columns: ", 
       paste(missing_cols, collapse = ", "), 
       "\nPlease ensure your TSV contains 'SampleID', 'Group', 'Day', and 'SubjectID'.")
}

# Align metadata and sequence table sample names
common_samples <- intersect(sample_names_seqtab, metadata$SampleID)
if (length(common_samples) == 0) {
  stop("ERROR: No matching samples found between DADA2 output and metadata.")
}

cat("Found", length(common_samples), "matching samples between ASV table and metadata.\n")

# Subset sequence table and metadata to match
seqtab <- seqtab[common_samples, , drop = FALSE]
metadata <- metadata %>% filter(SampleID %in% common_samples) %>% as.data.frame()
rownames(metadata) <- metadata$SampleID

# -------------------------------------------------------------------------
# Step 2: Build Phyloseq Object
# -------------------------------------------------------------------------
cat(">>> Constructing Phyloseq object...\n")
OTU <- otu_table(seqtab, taxa_are_rows = FALSE)
TAX <- tax_table(taxa)
META <- sample_data(metadata)

# Create raw phyloseq object
physeq <- phyloseq(OTU, TAX, META)

# Standardize tax table names (replace NA or empty assignments with placeholder)
tax_table(physeq)[is.na(tax_table(physeq))] <- "Unclassified"
for (i in 1:ncol(tax_table(physeq))) {
  tax_table(physeq)[tax_table(physeq)[, i] == "", i] <- "Unclassified"
}

# Save phyloseq RDS object
saveRDS(physeq, file.path(PHYLO_DIR, "phyloseq_obj.rds"))
cat("Phyloseq object saved to: ", file.path(PHYLO_DIR, "phyloseq_obj.rds"), "\n")

# -------------------------------------------------------------------------
# Step 3: Export Cleaned Taxonomic Abundance Tables
# -------------------------------------------------------------------------
cat(">>> Exporting abundance tables at different taxonomic ranks...\n")

wb <- createWorkbook()

# We will export tables for: ASV, Phylum, Class, Order, Family, Genus
ranks <- c("ASV", "Phylum", "Class", "Order", "Family", "Genus")

for (rank in ranks) {
  cat("Processing rank: ", rank, "\n")
  
  if (rank == "ASV") {
    # Extract ASV table and attach taxonomy info
    counts <- as.data.frame(t(otu_table(physeq)))
    tax_info <- as.data.frame(tax_table(physeq))
    merged_tab <- cbind(ASV_ID = rownames(counts), tax_info, counts)
  } else {
    # Collapse (glom) by specified rank
    glom_physeq <- tax_glom(physeq, taxrank = rank, NArm = FALSE)
    counts <- as.data.frame(t(otu_table(glom_physeq)))
    tax_info <- as.data.frame(tax_table(glom_physeq))[, 1:match(rank, colnames(tax_table(physeq))), drop = FALSE]
    
    # Merge taxonomy and counts
    merged_tab <- cbind(tax_info, counts)
    rownames(merged_tab) <- NULL
    
    # Group by taxonomy and sum to resolve duplicate names
    merged_tab <- merged_tab %>%
      group_by(across(1:all_of(rank))) %>%
      summarise(across(everything(), sum), .groups = "drop")
  }
  
  # Save as TSV
  write_tsv(merged_tab, file.path(PHYLO_DIR, paste0("abundance_", tolower(rank), ".tsv")))
  
  # Add to Excel workbook (sheets must be <= 31 chars)
  addWorksheet(wb, rank)
  writeData(wb, rank, merged_tab)
}

# Save Excel Workbook
saveWorkbook(wb, file.path(PHYLO_DIR, "taxa_abundance_tables.xlsx"), overwrite = TRUE)

# Export a dedicated pure ASV count table for PICRUSt2
cat(">>> Exporting standard ASV count table for PICRUSt2...\n")
asv_counts_picrust <- as.data.frame(t(otu_table(physeq)))
asv_counts_picrust <- cbind("#OTU ID" = rownames(asv_counts_picrust), asv_counts_picrust)
write_tsv(asv_counts_picrust, file.path(PHYLO_DIR, "asv_table_picrust.tsv"))

cat("=========================================================================\n")
cat("Phyloseq Building & Export Completed successfully!\n")
cat("Excel Workbook: ", file.path(PHYLO_DIR, "taxa_abundance_tables.xlsx"), "\n")
cat("Phyloseq RDS:   ", file.path(PHYLO_DIR, "phyloseq_obj.rds"), "\n")
cat("=========================================================================\n")
