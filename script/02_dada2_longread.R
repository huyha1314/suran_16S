#!/usr/bin/env Rscript

# =========================================================================
#       Long-read 16S Taxonomic Profiling Pipeline (Emu & Kraken2/Bracken)
# =========================================================================

cat("Loading required libraries...\n")
library(tidyverse)

# Retrieve execution parameters
THREADS <- as.numeric(Sys.getenv("THREADS", "64"))
TRIMMED_DIR <- Sys.getenv("TRIMMED_DIR", "./results/01_trimmed")
DADA2_DIR <- Sys.getenv("DADA2_DIR", "./results/02_dada2")
DB_DIR <- Sys.getenv("DB_DIR", "./data/db")

EMU_DB <- file.path(DB_DIR, "emu")
KRAKEN_DB <- file.path(DB_DIR, "kraken2/16S_SILVA138_k2db")

cat("=========================================================================\n")
cat("Parameters:\n")
cat("Trimmed Fastq Dir:       ", TRIMMED_DIR, "\n")
cat("Output Dir:              ", DADA2_DIR, "\n")
cat("Emu Database Dir:        ", EMU_DB, "\n")
cat("Kraken2 Database Dir:    ", KRAKEN_DB, "\n")
cat("Threads:                 ", THREADS, "\n")
cat("=========================================================================\n")

# Verify databases exist
if (!dir.exists(EMU_DB) || !file.exists(file.path(EMU_DB, "species_taxid.fasta"))) {
  stop("ERROR: Emu reference database not found at: ", EMU_DB)
}
if (!dir.exists(KRAKEN_DB) || !file.exists(file.path(KRAKEN_DB, "hash.k2d"))) {
  stop("ERROR: Kraken2 reference database not found at: ", KRAKEN_DB)
}

# Create output subdirectories
emu_out_dir <- file.path(DADA2_DIR, "emu")
kraken_out_dir <- file.path(DADA2_DIR, "kraken2")
dir.create(emu_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(kraken_out_dir, showWarnings = FALSE, recursive = TRUE)

# Identify verified trimmed fastq files
fastq_files <- list.files(TRIMMED_DIR, pattern = "_trimmed\\.fastq\\.gz$", full.names = FALSE)
if (length(fastq_files) == 0) {
  stop("ERROR: No trimmed FASTQ files (*_trimmed.fastq.gz) found in: ", TRIMMED_DIR)
}

sample_ids <- str_remove(fastq_files, "_trimmed\\.fastq\\.gz$")
cat("Processing", length(sample_ids), "samples:\n")
print(sample_ids)

# -------------------------------------------------------------------------
# Step 1: Run Emu taxonomic profiling sequentially
# -------------------------------------------------------------------------
cat("\n>>> Running Emu taxonomic classification on all samples...\n")
for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  fastq_path <- file.path(TRIMMED_DIR, fastq_files[i])
  rel_abund_file <- file.path(emu_out_dir, paste0(sample_id, "_rel-abundance.tsv"))
  
  cat("\n⏳ [Emu] Processing sample:", sample_id, "(", i, "/", length(sample_ids), ")\n")
  
  if (file.exists(rel_abund_file)) {
    cat("  [SKIP] Emu taxonomic profile already exists at:", rel_abund_file, "\n")
  } else {
    # Run Emu abundance
    cmd <- sprintf("pixi run emu abundance %s --db %s --threads %d --output-dir %s --output-basename %s --keep-counts --type map-ont",
                   fastq_path, EMU_DB, THREADS, emu_out_dir, sample_id)
    cat("Running:", cmd, "\n")
    system(cmd)
  }
}

# Combine Emu outputs
cat("\n>>> Combining Emu abundance estimates at species level...\n")
combine_cmd <- sprintf("pixi run emu combine-outputs --split-tables --counts %s species", emu_out_dir)
cat("Running:", combine_cmd, "\n")
system(combine_cmd)

# -------------------------------------------------------------------------
# Step 2: Run Kraken2 and Bracken sequentially
# -------------------------------------------------------------------------
cat("\n>>> Running Kraken2 & Bracken classification on all samples...\n")
bracken_files <- c()
bracken_names <- c()

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  fastq_path <- file.path(TRIMMED_DIR, fastq_files[i])
  
  k2_report <- file.path(kraken_out_dir, paste0(sample_id, "_kraken2.report"))
  k2_out <- file.path(kraken_out_dir, paste0(sample_id, "_kraken2.out"))
  bracken_out <- file.path(kraken_out_dir, paste0(sample_id, "_bracken.out"))
  bracken_report <- file.path(kraken_out_dir, paste0(sample_id, "_bracken.report"))
  
  cat("\n⏳ [Kraken2/Bracken] Processing sample:", sample_id, "(", i, "/", length(sample_ids), ")\n")
  
  if (file.exists(bracken_out)) {
    cat("  [SKIP] Kraken2/Bracken taxonomic profile already exists at:", bracken_out, "\n")
  } else {
    # Run Kraken2
    k2_cmd <- sprintf("pixi run kraken2 --db %s --threads %d --gzip-compressed %s --report %s --output %s",
                      KRAKEN_DB, THREADS, fastq_path, k2_report, k2_out)
    cat("Running:", k2_cmd, "\n")
    system(k2_cmd)
    
    # Run Bracken (genus level, standard read length 150)
    bracken_cmd <- sprintf("pixi run bracken -d %s -i %s -o %s -w %s -r 150 -l G",
                           KRAKEN_DB, k2_report, bracken_out, bracken_report)
    cat("Running:", bracken_cmd, "\n")
    system(bracken_cmd)
  }
  
  # Track for combine script
  if (file.exists(bracken_out)) {
    bracken_files <- c(bracken_files, bracken_out)
    bracken_names <- c(bracken_names, sample_id)
  }
}

# Combine Bracken outputs if files exist
combined_bracken_file <- file.path(kraken_out_dir, "bracken_combined.tsv")
if (length(bracken_files) > 0) {
  cat("\n>>> Combining Bracken abundance estimates...\n")
  combine_bracken_cmd <- sprintf("pixi run combine_bracken_outputs.py --files %s --names %s -o %s",
                                 paste(bracken_files, collapse = " "),
                                 paste(bracken_names, collapse = ","),
                                 combined_bracken_file)
  cat("Running:", combine_bracken_cmd, "\n")
  system(combine_bracken_cmd)
}

# -------------------------------------------------------------------------
# Step 3: Parse Emu outputs and build standard DADA2-style RDS objects
# -------------------------------------------------------------------------
cat("\n>>> Compiling primary Emu taxonomy and abundance RDS files...\n")
emu_abundance_file <- file.path(emu_out_dir, "emu-combined-abundance-species-counts.tsv")
emu_taxonomy_file <- file.path(emu_out_dir, "emu-combined-taxonomy-species.tsv")

if (!file.exists(emu_abundance_file) || !file.exists(emu_taxonomy_file)) {
  stop("ERROR: Emu combined output files were not generated!")
}

# Read Emu abundance matrix (taxa x samples)
emu_counts <- read.delim(emu_abundance_file, sep="\t", header=TRUE, check.names=FALSE, row.names=1)
emu_counts <- emu_counts[rownames(emu_counts) != "" & !is.na(rownames(emu_counts)), , drop = FALSE]

# Read Emu taxonomy lineage
emu_tax_df <- read.delim(emu_taxonomy_file, sep="\t", header=TRUE, check.names=FALSE, row.names=1)
emu_tax_df <- emu_tax_df[rownames(emu_tax_df) != "" & !is.na(rownames(emu_tax_df)), , drop = FALSE]

# Format to standard DADA2 structures:
# seqtab: samples x ASVs (where ASV is the tax_id)
seqtab_emu <- t(as.matrix(emu_counts))

# taxonomy: ASVs x 7 ranks (Kingdom, Phylum, Class, Order, Family, Genus, Species)
taxa_emu <- emu_tax_df %>%
  select(superkingdom, phylum, class, order, family, genus) %>%
  rename(Kingdom = superkingdom, Phylum = phylum, Class = class, Order = order,
         Family = family, Genus = genus) %>%
  mutate(Species = rownames(emu_tax_df)) %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  as.matrix()
rownames(taxa_emu) <- rownames(emu_tax_df)

# Align row/column names
common_taxa <- intersect(colnames(seqtab_emu), rownames(taxa_emu))
common_taxa <- common_taxa[common_taxa != "" & !is.na(common_taxa)]
seqtab_emu <- seqtab_emu[, common_taxa, drop = FALSE]
taxa_emu <- taxa_emu[common_taxa, , drop = FALSE]

# Save primary DADA2-style files to the main results dir so downstream pipelines consume them directly!
saveRDS(seqtab_emu, file.path(DADA2_DIR, "seqtab_nochim.rds"))
saveRDS(taxa_emu, file.path(DADA2_DIR, "taxonomy.rds"))

# -------------------------------------------------------------------------
# Step 4: Parse Bracken outputs and build standard RDS files
# -------------------------------------------------------------------------
if (file.exists(combined_bracken_file)) {
  cat("\n>>> Compiling alternative Bracken taxonomy and abundance RDS files...\n")
  bracken_df <- read.delim(combined_bracken_file, sep="\t", header=TRUE, check.names=FALSE)
  
  # Bracken columns: name, taxonomy_id, taxonomy_lvl, then samples
  bracken_counts <- bracken_df[, 4:ncol(bracken_df), drop = FALSE]
  rownames(bracken_counts) <- bracken_df$taxonomy_id
  seqtab_bracken <- t(as.matrix(bracken_counts))
  
  # Map taxonomy for Bracken using Emu lineage details or defaults
  taxa_bracken <- matrix("Unclassified", nrow = nrow(bracken_df), ncol = 7)
  rownames(taxa_bracken) <- bracken_df$taxonomy_id
  colnames(taxa_bracken) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # Populate matched species/genera from Bracken names
  for (j in 1:nrow(bracken_df)) {
    taxid_str <- as.character(bracken_df$taxonomy_id[j])
    spec_name <- bracken_df$name[j]
    
    # If present in Emu, copy full lineage
    if (taxid_str %in% rownames(taxa_emu)) {
      taxa_bracken[taxid_str, ] <- taxa_emu[taxid_str, ]
    } else {
      # Fallback parsing for Genus-level Bracken
      taxa_bracken[taxid_str, "Kingdom"] <- "Bacteria"
      taxa_bracken[taxid_str, "Genus"] <- spec_name
      taxa_bracken[taxid_str, "Species"] <- "sp."
    }
  }
  
  saveRDS(seqtab_bracken, file.path(kraken_out_dir, "seqtab_nochim.rds"))
  saveRDS(taxa_bracken, file.path(kraken_out_dir, "taxonomy.rds"))
}

# -------------------------------------------------------------------------
# Step 5: Read Tracking Summary
# -------------------------------------------------------------------------
cat("\n>>> Compiling read-tracking summary...\n")
input_reads <- numeric(length(sample_ids))
names(input_reads) <- sample_ids

# Fastq reading count using basic zcat/gzip linecount
for (i in seq_along(sample_ids)) {
  fastq_path <- file.path(TRIMMED_DIR, fastq_files[i])
  linecount <- as.numeric(system(sprintf("zcat %s | wc -l", fastq_path), intern = TRUE))
  input_reads[sample_ids[i]] <- floor(linecount / 4)
}

# Sum of species-assigned reads from Emu counts
emu_assigned_reads <- rowSums(seqtab_emu)

# Create simplified read tracking table
track <- data.frame(
  input = input_reads,
  filtered = input_reads,
  denoised = emu_assigned_reads,
  non_chimera = emu_assigned_reads,
  row.names = sample_ids
)

write.table(track, file.path(DADA2_DIR, "read_tracking.tsv"), 
            sep = "\t", quote = FALSE, col.names = NA)

cat("\nRead tracking table:\n")
print(track)

cat("\n=========================================================================\n")
cat("Long-Read Emu & Kraken2/Bracken Profiling Completed Successfully!\n")
cat("Emu Counts Matrix:   ", file.path(DADA2_DIR, "seqtab_nochim.rds"), "\n")
cat("Emu Taxonomy Table:  ", file.path(DADA2_DIR, "taxonomy.rds"), "\n")
cat("Kraken2/Bracken Dir: ", kraken_out_dir, "\n")
cat("Read Tracking TSV:   ", file.path(DADA2_DIR, "read_tracking.tsv"), "\n")
cat("=========================================================================\n")
