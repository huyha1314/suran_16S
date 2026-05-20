#!/usr/bin/env Rscript

# =========================================================================
#                 16S rRNA V3-V4 DADA2 Denoising Pipeline
# =========================================================================

cat("Loading required libraries...\n")
library(dada2)
library(ggplot2)

# Get parameters from environment variables (with default fallbacks)
TRIM_DIR <- Sys.getenv("TRIM_DIR", "./results/01_trimmed")
DADA2_DIR <- Sys.getenv("DADA2_DIR", "./results/02_dada2")
DB_DIR <- Sys.getenv("DB_DIR", "./data/db")
THREADS <- as.integer(Sys.getenv("THREADS", "8"))

# For V3-V4 300 PE, default truncations and maxEE
# Sum of truncLen (270 + 210 = 480 bp) must be greater than target trimmed length (~420 bp) + 20 bp overlap.
# 480 bp is plenty of overlap (~60 bp) even with low-quality reverse ends trimmed at 210.
TRUNC_LEN_F <- as.integer(Sys.getenv("TRUNC_LEN_F", "270"))
TRUNC_LEN_R <- as.integer(Sys.getenv("TRUNC_LEN_R", "210"))
MAX_EE_F <- as.numeric(Sys.getenv("MAX_EE_F", "2"))
MAX_EE_R <- as.numeric(Sys.getenv("MAX_EE_R", "5")) # Relaxed for PE300 reverse reads

cat("=========================================================================\n")
cat("Parameters:\n")
cat("Trimmed Fastq Dir:      ", TRIM_DIR, "\n")
cat("DADA2 Output Dir:       ", DADA2_DIR, "\n")
cat("Database Dir:           ", DB_DIR, "\n")
cat("Threads:                ", THREADS, "\n")
cat("Truncation Length F/R:  ", TRUNC_LEN_F, " / ", TRUNC_LEN_R, "\n")
cat("Max Expected Errors F/R:", MAX_EE_F, " / ", MAX_EE_R, "\n")
cat("=========================================================================\n")

# Create output directories
dir.create(DADA2_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DB_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Step 0: Verify Reference Databases
# -------------------------------------------------------------------------
silva_train_file <- file.path(DB_DIR, "silva_nr99_v138.1_train_set.fa.gz")
silva_species_file <- file.path(DB_DIR, "silva_species_assignment_v138.1.fa.gz")

if (!file.exists(silva_train_file) || !file.exists(silva_species_file)) {
  stop("ERROR: SILVA reference database files are missing in ", DB_DIR, 
       "\nPlease download them first by running: pixi run download_db")
}
cat("Database verification complete: SILVA training set and species database verified.\n")

# -------------------------------------------------------------------------
# Step 1: Read Fastq Files from Sample TSV
# -------------------------------------------------------------------------
sample_tsv_path <- Sys.getenv("SAMPLE_TSV", "./sample.tsv")
if (!file.exists(sample_tsv_path)) {
  stop("ERROR: Sample sheet not found at: ", sample_tsv_path, "\n",
       "Please generate it first by running: pixi run generate_samples")
}

cat("Reading sample sheet from:", sample_tsv_path, "...\n")
df_samples <- read.delim(sample_tsv_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (!"SampleID" %in% colnames(df_samples)) {
  stop("ERROR: The sample sheet must contain a 'SampleID' column.")
}

sample.names <- df_samples$SampleID

# Construct paired trimmed FASTQ paths
fnFs <- file.path(TRIM_DIR, paste0(sample.names, "_trimmed_R1.fastq.gz"))
fnRs <- file.path(TRIM_DIR, paste0(sample.names, "_trimmed_R2.fastq.gz"))

# Verify trimmed files exist
existing_indices <- file.exists(fnFs) & file.exists(fnRs)
if (sum(existing_indices) == 0) {
  stop("ERROR: No trimmed FASTQ files found at: ", TRIM_DIR, 
       "\nPlease run Stage 1 (QC & Primer Trimming) first.")
}

# Keep only samples with existing trimmed FASTQs
sample.names <- sample.names[existing_indices]
fnFs <- fnFs[existing_indices]
fnRs <- fnRs[existing_indices]

cat("Processing", length(sample.names), "samples with verified trimmed FASTQs:\n")
print(sample.names)

# -------------------------------------------------------------------------
# Step 2: Quality Profiling (Before Denoising)
# -------------------------------------------------------------------------
cat(">>> Generating Quality Profile Plots...\n")
p_f <- plotQualityProfile(fnFs[1:min(length(fnFs), 2)])
p_r <- plotQualityProfile(fnRs[1:min(length(fnRs), 2)])
ggsave(file.path(DADA2_DIR, "quality_profile_R1_raw.png"), plot = p_f, width = 8, height = 6)
ggsave(file.path(DADA2_DIR, "quality_profile_R2_raw.png"), plot = p_r, width = 8, height = 6)

# -------------------------------------------------------------------------
# Step 3: Filter and Trim
# -------------------------------------------------------------------------
cat(">>> Filtering and trimming reads...\n")
filtFs <- file.path(DADA2_DIR, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(DADA2_DIR, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(TRUNC_LEN_F, TRUNC_LEN_R),
                     maxN = 0, 
                     maxEE = c(MAX_EE_F, MAX_EE_R), 
                     truncQ = 2, 
                     rm.phix = TRUE,
                     compress = TRUE, 
                     multithread = THREADS,
                     verbose = TRUE)

cat("Filtering stats summary:\n")
print(out)

# -------------------------------------------------------------------------
# Step 4: Learn Error Rates
# -------------------------------------------------------------------------
cat(">>> Learning error rates...\n")
errF <- learnErrors(filtFs, multithread = THREADS, verbose = TRUE)
errR <- learnErrors(filtRs, multithread = THREADS, verbose = TRUE)

# Save error rate plots
pdf(file.path(DADA2_DIR, "error_rates_F.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(file.path(DADA2_DIR, "error_rates_R.pdf"))
plotErrors(errR, nominalQ = TRUE)
dev.off()

# -------------------------------------------------------------------------
# Step 5: Denoising (ASV Inference)
# -------------------------------------------------------------------------
cat(">>> Running DADA2 denoising algorithm (ASV inference)...\n")
dadaFs <- dada(filtFs, err = errF, multithread = THREADS, pool = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = THREADS, pool = FALSE)

# -------------------------------------------------------------------------
# Step 6: Merge Paired-End Reads
# -------------------------------------------------------------------------
cat(">>> Merging read pairs...\n")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, 
                      minOverlap = 20, 
                      maxMismatch = 0, 
                      verbose = TRUE)

# -------------------------------------------------------------------------
# Step 7: Construct Sequence Table (ASVs)
# -------------------------------------------------------------------------
cat(">>> Constructing ASV table...\n")
seqtab <- makeSequenceTable(mergers)
cat("Initial ASV table dimension (samples x ASVs):", dim(seqtab), "\n")

# Distribution of sequence lengths
cat("ASV sequence length distribution:\n")
print(table(nchar(getSequences(seqtab))))

# -------------------------------------------------------------------------
# Step 8: Remove Chimeras
# -------------------------------------------------------------------------
cat(">>> Removing chimeras (Bimeras)...\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "consensus", 
                                    multithread = THREADS, 
                                    verbose = TRUE)
cat("Dimension after chimera removal (samples x ASVs):", dim(seqtab.nochim), "\n")
cat("Percentage of non-chimeric reads kept: ", 
    round(100 * sum(seqtab.nochim) / sum(seqtab), 2), "%\n")

# -------------------------------------------------------------------------
# Step 9: Taxonomic Classification (SILVA v138.1)
# -------------------------------------------------------------------------
cat(">>> Assigning taxonomy using SILVA training set...\n")
taxa <- assignTaxonomy(seqtab.nochim, silva_train_file, multithread = THREADS)

cat(">>> Adding species-level assignment in parallel...\n")
taxa <- addSpecies(taxa, silva_species_file, multithread = THREADS)

# -------------------------------------------------------------------------
# Step 9.5: Standardize ASV Names (Crucial for BIOM & PICRUSt2 Compatibility)
# -------------------------------------------------------------------------
cat(">>> Renaming ASV sequences to standard clean identifiers (ASV_00001 style)...\n")
asv_seqs <- colnames(seqtab.nochim)
asv_names <- sprintf("ASV_%05d", seq_along(asv_seqs))

# Rename count table columns and taxonomy table rows
colnames(seqtab.nochim) <- asv_names
rownames(taxa) <- asv_names

# Save Clean tables to RDS
saveRDS(seqtab.nochim, file.path(DADA2_DIR, "seqtab_nochim.rds"))
saveRDS(taxa, file.path(DADA2_DIR, "taxonomy.rds"))

# Export standard ASV Fasta matching the table IDs
asv_headers <- paste0(">", asv_names)
asv_fasta <- c(rbind(asv_headers, asv_seqs))
writeLines(asv_fasta, file.path(DADA2_DIR, "asvs.fasta"))

# -------------------------------------------------------------------------
# Step 10: Track Reads through Pipeline
# -------------------------------------------------------------------------
cat(">>> Compiling read-tracking summary...\n")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoised_F", "denoised_R", "merged", "non_chimera")
rownames(track) <- sample.names
write.table(track, file.path(DADA2_DIR, "read_tracking.tsv"), 
            sep = "\t", quote = FALSE, col.names = NA)

cat("Read tracking table:\n")
print(track)

cat("=========================================================================\n")
cat("DADA2 Denoising and Taxonomic Classification Completed Successfully!\n")
cat("ASV Table RDS:     ", file.path(DADA2_DIR, "seqtab_nochim.rds"), "\n")
cat("Taxonomy RDS:      ", file.path(DADA2_DIR, "taxonomy.rds"), "\n")
cat("ASV FASTA:         ", file.path(DADA2_DIR, "asvs.fasta"), "\n")
cat("Read Tracking TSV: ", file.path(DADA2_DIR, "read_tracking.tsv"), "\n")
cat("=========================================================================\n")
