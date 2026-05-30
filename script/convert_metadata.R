#!/usr/bin/env Rscript

cat("Loading R packages for longitudinal metadata generation...\n")
library(dplyr)
library(stringr)
library(readr)

# Read environment variables or use default paths
sample_tsv_path <- Sys.getenv("SAMPLE_TSV", "./sample.tsv")
metadata_tsv_path <- Sys.getenv("METADATA_PATH", "./metadata.tsv")

if (!file.exists(sample_tsv_path)) {
  stop("ERROR: Sample sheet not found at: ", sample_tsv_path, "\n",
       "Please generate it first by running: pixi run generate_samples")
}

cat("Reading sample sheet from:", sample_tsv_path, "...\n")
df_samples <- read_tsv(sample_tsv_path)

if (!"SampleID" %in% colnames(df_samples)) {
  stop("ERROR: The sample sheet must contain a 'SampleID' column.")
}

cat("Inferring longitudinal study groups...\n")
parsed_metadata <- df_samples %>%
  mutate(
    # 1. Dynamic prefix extraction with multi-tier fallback:
    # Tier 1: Extract everything before the last dash/underscore + digits (e.g. EM_10 -> EM, DC3_01 -> DC3, 0-12 -> 0)
    GroupPrefix = str_match(SampleID, "^(.+?)[_-]\\d+$")[,2],
    
    # Tier 2: Extract everything before the first dash/underscore (e.g. DP4-B9V2 -> DP4) if Tier 1 returned NA
    GroupPrefix = ifelse(is.na(GroupPrefix), str_match(SampleID, "^([A-Za-z0-9]+?)[_-]")[,2], GroupPrefix),
    
    # Tier 3: Default to full SampleID if no separator is found
    GroupPrefix = ifelse(is.na(GroupPrefix), SampleID, GroupPrefix),
    
    # 2. Standardize Groups dynamically (maintaining explicit groups for shortread baseline datasets)
    Group = case_when(
      GroupPrefix == "0"  ~ "Group1",
      GroupPrefix == "21" ~ "Group2",
      GroupPrefix == "35" ~ "Group3",
      is.na(GroupPrefix)  ~ "Group1",
      TRUE ~ GroupPrefix
    ),
    
    # Add sequential sample counter (STT)
    STT = row_number()
  ) %>%
  select(SampleID, Group, STT)

# Ensure data folder exists
dir.create(dirname(metadata_tsv_path), showWarnings = FALSE, recursive = TRUE)

# Save standardized TSV
write_tsv(parsed_metadata, metadata_tsv_path)

cat("SUCCESS: Automated metadata generated and saved to:", metadata_tsv_path, "\n")
print(head(parsed_metadata, 10))

