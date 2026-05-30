#!/usr/bin/env Rscript

library(phyloseq)
library(decontam)
library(ggplot2)

PHYLO_DIR <- Sys.getenv("PHYLO_DIR", "./results/03_phyloseq")
ps_path <- file.path(PHYLO_DIR, "phyloseq_obj.rds")

if (!file.exists(ps_path)) {
  stop("Phyloseq object not found. Run 03_phyloseq_prep.R first.")
}

ps <- readRDS(ps_path)
cat(">>> Starting Decontam Analysis...\n")

# Check if a control column exists in metadata
meta <- as(sample_data(ps), "data.frame")
control_col <- "is_control"

if (!control_col %in% colnames(meta)) {
  cat("⚠️ WARNING: No 'is_control' column found in metadata.tsv. Skipping Decontam.\n")
  cat("To use decontam, add an 'is_control' column (TRUE for blanks/negatives, FALSE for samples).\n")
  quit(save = "no", status = 0)
}

# Run Decontam using Prevalence method (comparing samples vs negative controls)
cat(">>> Identifying contaminants using negative controls (prevalence method)...\n")
sample_data(ps)$is_control <- as.logical(sample_data(ps)[[control_col]])
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is_control", threshold=0.1)

table(contamdf.prev$contaminant)
contaminant_taxa <- rownames(contamdf.prev[contamdf.prev$contaminant, ])
cat("Found", length(contaminant_taxa), "contaminant taxa.\n")

# Filter contaminants out
ps_clean <- prune_taxa(!taxa_names(ps) %in% contaminant_taxa, ps)

# Save cleaned phyloseq object
saveRDS(ps_clean, file.path(PHYLO_DIR, "phyloseq_obj.rds"))
cat(">>> Decontam complete! Cleaned phyloseq object saved.\n")
