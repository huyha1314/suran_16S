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

cat("Inferring longitudinal study groups, days, and SubjectIDs...\n")
parsed_metadata <- df_samples %>%
  mutate(
    # Extract the timepoint suffix (e.g. 21-3.2 -> Day2)
    DayNum = str_match(SampleID, "\\.(\\d+)$")[,2],
    Day = ifelse(!is.na(DayNum), paste0("Day", DayNum), "Day1"),
    
    # Extract prefix subject codes (e.g. 21-3.2 -> Subject_21_3)
    SubjPrefix = str_match(SampleID, "^(\\d+)-(\\d+)\\.")[,2],
    SubjNum = str_match(SampleID, "^(\\d+)-(\\d+)\\.")[,3],
    SubjectID = ifelse(!is.na(SubjPrefix) & !is.na(SubjNum), 
                       paste0("Subject_", SubjPrefix, "_", SubjNum), 
                       paste0("Subject_", SampleID)),
    
    # Automatically infer experimental groups based on the subject prefix range:
    # Prefix '0' -> Group1 (Nhóm 1)
    # Prefix '21' -> Group2 (Nhóm 2)
    # Prefix '35' -> Group3 (Nhóm 3)
    Group = case_when(
      SubjPrefix == "0"  ~ "Group1",
      SubjPrefix == "21" ~ "Group2",
      SubjPrefix == "35" ~ "Group3",
      TRUE ~ ifelse(!is.na(SubjPrefix), paste0("Group_", SubjPrefix), "Group1")
    ),
    
    # Add sequential sample counter (STT)
    STT = row_number()
  ) %>%
  select(SampleID, Group, Day, SubjectID, STT)

# Ensure data folder exists
dir.create(dirname(metadata_tsv_path), showWarnings = FALSE, recursive = TRUE)

# Save standardized TSV
write_tsv(parsed_metadata, metadata_tsv_path)

cat("SUCCESS: Automated metadata generated and saved to:", metadata_tsv_path, "\n")
print(head(parsed_metadata, 10))
