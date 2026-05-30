#!/usr/bin/env Rscript

# =========================================================================
#                 16S rRNA Phyloseq Object Generation & Export
# =========================================================================

cat("Loading required libraries...\n")
library(phyloseq)
library(tidyverse)
library(openxlsx)
library(DECIPHER)
library(Biostrings)
library(ape)

# Get parameters from environment variables (with default fallbacks)
DADA2_DIR <- Sys.getenv("DADA2_DIR", "./results/02_dada2")
PHYLO_DIR <- Sys.getenv("PHYLO_DIR", "./results/03_phyloseq")
METADATA_PATH <- Sys.getenv("METADATA_PATH", "./metadata.tsv")
MODE <- Sys.getenv("MODE", "shortread")

cat("=========================================================================\n")
cat("Parameters:\n")
cat("DADA2 Output Dir:   ", DADA2_DIR, "\n")
cat("Phyloseq Output Dir:", PHYLO_DIR, "\n")
cat("Metadata File Path: ", METADATA_PATH, "\n")
cat("Execution Mode:     ", MODE, "\n")
cat("=========================================================================\n")

# Create output directory
dir.create(PHYLO_DIR, showWarnings = FALSE, recursive = TRUE)
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", "./results")

if (MODE != "longread") {
  # Load DADA2 outputs
  seqtab_path <- file.path(DADA2_DIR, "seqtab_nochim.rds")
  taxa_path <- file.path(DADA2_DIR, "taxonomy.rds")
  
  if (!file.exists(seqtab_path) || !file.exists(taxa_path)) {
    stop("ERROR: DADA2 output files not found. Please run DADA2 step first.")
  }
  
  seqtab <- readRDS(seqtab_path)
  taxa <- readRDS(taxa_path)
  
  sample_names_seqtab <- rownames(seqtab)
}

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
required_cols <- c("SampleID", "Group")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("ERROR: Standardized metadata is missing required columns: ", 
       paste(missing_cols, collapse = ", "), 
       "\nPlease ensure your TSV contains 'SampleID' and 'Group'.")
}

if (MODE != "longread") {
  # Align metadata and sequence table sample names
  common_samples <- intersect(sample_names_seqtab, metadata$SampleID)
  if (length(common_samples) == 0) {
    stop("ERROR: No matching samples found between DADA2 output and metadata.")
  }
  
  cat("Found", length(common_samples), "matching samples between ASV table and metadata.\n")
  
  # Subset sequence table and metadata to match
  seqtab <- seqtab[common_samples, , drop = FALSE]
  metadata <- metadata %>% filter(SampleID %in% common_samples) %>% as.data.frame()
} else {
  metadata <- as.data.frame(metadata)
}
rownames(metadata) <- metadata$SampleID
sample_data_obj <- sample_data(metadata)

# -------------------------------------------------------------------------
# Step 2: Build Phylogenetic Tree & Phyloseq Object
# -------------------------------------------------------------------------
if (MODE == "longread") {
  cat(">>> Long-read mode detected. Bypassing DECIPHER/FastTree & DADA2 entirely.\n")
  
  # =========================================================================
  # 1. EMU PIPELINE TO PHYLOSEQ
  # =========================================================================
  cat(">>> Loading Emu Data...\n")
  
  emu_abundance_file <- file.path(RESULTS_DIR, "emu_combined", "emu-combined-abundance.tsv")
  if (!file.exists(emu_abundance_file)) {
    stop("ERROR: Emu combined abundance file not found at: ", emu_abundance_file)
  }
  
  emu_data <- read_tsv(emu_abundance_file)
  
  emu_tax <- emu_data %>% dplyr::select(superkingdom, phylum, class, order, family, genus, species) %>%
    dplyr::rename(Kingdom = superkingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, "Unclassified")))
  emu_tax_mat <- as.matrix(emu_tax)
  rownames(emu_tax_mat) <- make.unique(emu_tax$Species)
  
  emu_counts <- emu_data %>% dplyr::select(-superkingdom, -phylum, -class, -order, -family, -genus, -species) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0)))
  emu_count_mat <- round(as.matrix(emu_counts))
  rownames(emu_count_mat) <- make.unique(emu_tax$Species)
  
  ps_emu <- phyloseq(otu_table(emu_count_mat, taxa_are_rows = TRUE),
                     tax_table(emu_tax_mat),
                     sample_data_obj)
  
  ps_emu_phylum <- tax_glom(ps_emu, taxrank = "Phylum")
  ps_emu_rel <- transform_sample_counts(ps_emu_phylum, function(x) x / sum(x))
  
  plot_emu <- plot_bar(ps_emu_rel, fill = "Phylum") + 
    geom_bar(stat = "identity") + 
    ggtitle("Taxonomic Profile (Emu)") + 
    theme_minimal()
  
  ggsave(file.path(PHYLO_DIR, "emu_barplot.png"), plot_emu, width = 10, height = 6)
  
  # =========================================================================
  # 2. BRACKEN PIPELINE TO PHYLOSEQ
  # =========================================================================
  cat(">>> Loading Bracken Data...\n")
  bracken_biom_file <- file.path(RESULTS_DIR, "bracken_combined.biom")
  
  if (file.exists(bracken_biom_file)) {
    ps_bracken <- import_biom(bracken_biom_file)
    sample_names(ps_bracken) <- stringr::str_remove(sample_names(ps_bracken), "_bracken$")
    sample_data(ps_bracken) <- sample_data_obj
    colnames(tax_table(ps_bracken)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    ps_bracken_phylum <- tax_glom(ps_bracken, taxrank = "Phylum")
    ps_bracken_rel <- transform_sample_counts(ps_bracken_phylum, function(x) x / sum(x))
    
    plot_bracken <- plot_bar(ps_bracken_rel, fill = "Phylum") + 
      geom_bar(stat = "identity") + 
      ggtitle("Taxonomic Profile (Bracken)") + 
      theme_minimal()
    
    ggsave(file.path(PHYLO_DIR, "bracken_barplot.png"), plot_bracken, width = 10, height = 6)
    
    # Save Bracken phyloseq object
    saveRDS(ps_bracken, file.path(PHYLO_DIR, "phyloseq_bracken.rds"))
  } else {
    cat("Warning: Bracken biom file not found, skipping Bracken plot.\n")
  }
  
  cat(">>> Success: Both Emu and Bracken barplots generated!\n")
  physeq <- ps_emu # Keep Emu as primary object downstream
} else {
  cat(">>> Aligning ASV sequences with DECIPHER...\n")
  asv_fasta_path <- file.path(DADA2_DIR, "asvs.fasta")
  seqs <- readDNAStringSet(asv_fasta_path)
  alignment <- AlignSeqs(seqs, anchor=NA, processors=1)
  
  cat(">>> Building phylogenetic tree with FastTree...\n")
  aligned_fasta <- file.path(PHYLO_DIR, "aligned_seqs.fasta")
  tree_out <- file.path(PHYLO_DIR, "tree.nwk")
  writeXStringSet(alignment, aligned_fasta)
  
  # Run FastTree (GTR+CAT model for nucleotides)
  exit_code <- system2("fasttree", args = c("-nt", "-gtr", aligned_fasta), stdout = tree_out)
  if (exit_code != 0) {
    stop("ERROR: FastTree failed to build the phylogenetic tree.")
  }
  tree <- read.tree(tree_out)
  
  cat(">>> Constructing Phyloseq object...\n")
  OTU <- otu_table(seqtab, taxa_are_rows = FALSE)
  TAX <- tax_table(taxa)
  META <- sample_data_obj
  
  # Create raw phyloseq object with tree
  physeq <- phyloseq(OTU, TAX, META, phy_tree(tree))
}

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

# We will export tables for: ASV (shortread only), Phylum, Class, Order, Family, Genus, Species
if (MODE == "longread") {
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
} else {
  ranks <- c("ASV", "Phylum", "Class", "Order", "Family", "Genus", "Species")
}

for (rank in ranks) {
  cat("Processing rank: ", rank, "\n")
  
  if (rank == "ASV") {
    # Extract ASV table and attach taxonomy info
    if (taxa_are_rows(physeq)) {
      counts <- as.data.frame(otu_table(physeq))
    } else {
      counts <- as.data.frame(t(otu_table(physeq)))
    }
    tax_info <- as.data.frame(tax_table(physeq))
    merged_tab <- cbind(ASV_ID = rownames(counts), tax_info, counts)
  } else {
    # Collapse (glom) by specified rank
    glom_physeq <- tax_glom(physeq, taxrank = rank, NArm = FALSE)
    if (taxa_are_rows(glom_physeq)) {
      counts <- as.data.frame(otu_table(glom_physeq))
    } else {
      counts <- as.data.frame(t(otu_table(glom_physeq)))
    }
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

# Export a dedicated pure ASV count table for PICRUSt2 (short-read only)
if (MODE != "longread") {
  cat(">>> Exporting standard ASV count table for PICRUSt2...\n")
  asv_counts_picrust <- as.data.frame(t(otu_table(physeq)))
  asv_counts_picrust <- cbind("#OTU ID" = rownames(asv_counts_picrust), asv_counts_picrust)
  write_tsv(asv_counts_picrust, file.path(PHYLO_DIR, "asv_table_picrust.tsv"))
}

cat("=========================================================================\n")
cat("Phyloseq Building & Export Completed successfully!\n")
cat("Excel Workbook: ", file.path(PHYLO_DIR, "taxa_abundance_tables.xlsx"), "\n")
cat("Phyloseq RDS:   ", file.path(PHYLO_DIR, "phyloseq_obj.rds"), "\n")
cat("=========================================================================\n")
