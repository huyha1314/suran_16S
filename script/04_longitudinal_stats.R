#!/usr/bin/env Rscript

# =========================================================================
#       16S rRNA Statistical Analysis & Visualization (Group-wise)
# =========================================================================

cat("Loading required libraries...\n")
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)

# Get parameters from environment variables (with default fallbacks)
PHYLO_DIR <- Sys.getenv("PHYLO_DIR", "./results/03_phyloseq")
STATS_DIR <- Sys.getenv("STATS_DIR", "./results/04_stats")

cat("=========================================================================\n")
cat("Parameters:\n")
cat("Phyloseq Dir: ", PHYLO_DIR, "\n")
cat("Stats Dir:    ", STATS_DIR, "\n")
cat("=========================================================================\n")

# Create output directories
dir.create(STATS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(STATS_DIR, "alpha"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(STATS_DIR, "beta"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(STATS_DIR, "differential"), showWarnings = FALSE, recursive = TRUE)

# Load Phyloseq object
physeq_path <- file.path(PHYLO_DIR, "phyloseq_obj.rds")
if (!file.exists(physeq_path)) {
  stop("ERROR: Phyloseq object not found at ", physeq_path)
}
physeq <- readRDS(physeq_path)

# Extract metadata
metadata <- as(sample_data(physeq), "data.frame")

# Validate metadata columns for group analysis
required_cols <- c("Group")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("ERROR: Metadata is missing required columns for analysis: ", 
       paste(missing_cols, collapse = ", "), 
       "\nEnsure 'Group' is populated.")
}

# Ensure factors are set correctly
metadata$Group <- factor(metadata$Group)

# -------------------------------------------------------------------------
# 1. Alpha Diversity Analysis (ANOVA & Kruskal-Wallis)
# -------------------------------------------------------------------------
cat(">>> Section 1: Running Alpha Diversity Group Comparison...\n")

cat(">>> Preparing OTU table for Alpha Diversity (Rounding to integers)...\n")
# 1. Force the OTU table to contain only whole integers
otu_table(physeq) <- round(otu_table(physeq))

# Compute alpha diversity metrics
alpha_meas <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Chao1"))
alpha_data <- cbind(metadata, alpha_meas)
write.table(alpha_data, file.path(STATS_DIR, "alpha", "alpha_diversity_table.tsv"), 
            sep = "\t", quote = FALSE, col.names = NA)

# Loop through Shannon and Observed metrics
metrics <- c("Shannon", "Observed")
for (metric in metrics) {
  cat("Comparing", metric, "across groups...\n")
  
  # Fit ANOVA
  formula_str <- paste(metric, "~ Group")
  fit_anova <- aov(as.formula(formula_str), data = alpha_data)
  anova_res <- summary(fit_anova)[[1]]
  write.table(as.data.frame(anova_res), 
              file.path(STATS_DIR, "alpha", paste0("anova_", tolower(metric), ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # Fit Kruskal-Wallis
  kw_res <- kruskal.test(as.formula(formula_str), data = alpha_data)
  kw_df <- data.frame(
    Statistic = kw_res$statistic,
    Parameter = kw_res$parameter,
    P_Value = kw_res$p.value,
    Method = kw_res$method
  )
  write.table(kw_df, 
              file.path(STATS_DIR, "alpha", paste0("kruskal_", tolower(metric), ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
}

# Rarefaction Curve Analysis omitted. Handled dynamically on the fly by vis/4.rarefaction.R


# -------------------------------------------------------------------------
# 2. Beta Diversity (PERMANOVA, PERMDISP, PCoA, NMDS)
# -------------------------------------------------------------------------
cat(">>> Section 2: Running Beta Diversity Analysis...\n")

physeq_prop <- transform_sample_counts(physeq, function(x) x / sum(x))

# Calculate distances (including UniFrac which requires the phylogenetic tree)
bray_dist <- phyloseq::distance(physeq_prop, method = "bray")
jaccard_dist <- phyloseq::distance(physeq_prop, method = "jaccard")

dists <- list("Bray-Curtis" = bray_dist, "Jaccard" = jaccard_dist)

# Attempt to calculate UniFrac if tree is present
if (!is.null(phy_tree(physeq_prop, errorIfNULL = FALSE))) {
  cat("Phylogenetic tree found. Calculating UniFrac distances...\n")
  wunifrac_dist <- phyloseq::distance(physeq_prop, method = "wunifrac")
  unifrac_dist <- phyloseq::distance(physeq_prop, method = "unifrac")
  dists[["Weighted-UniFrac"]] <- wunifrac_dist
  dists[["Unweighted-UniFrac"]] <- unifrac_dist
} else {
  cat("WARNING: No phylogenetic tree found. Skipping UniFrac distances.\n")
}

for (dist_name in names(dists)) {
  curr_dist <- dists[[dist_name]]
  clean_name <- tolower(gsub("-", "_", dist_name))
  
  cat("Analyzing", dist_name, "beta diversity...\n")
  
  # 1. Run PERMANOVA
  permanova_res <- adonis2(curr_dist ~ Group, 
                           data = metadata, 
                           permutations = 999)
  
  write.table(as.data.frame(permanova_res), 
              file.path(STATS_DIR, "beta", paste0("permanova_", clean_name, ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
              
  # 2. Run PERMDISP (betadisper)
  dispersion <- betadisper(curr_dist, metadata$Group)
  permdisp_res <- permutest(dispersion, permutations = 999)
  
  write.table(as.data.frame(permdisp_res$tab), 
              file.path(STATS_DIR, "beta", paste0("permdisp_", clean_name, ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # 3. Perform PCoA Ordination
  pcoa_ord <- ordinate(physeq_prop, method = "MDS", distance = curr_dist)
  
  pcoa_df <- data.frame(pcoa_ord$vectors[, 1:2])
  colnames(pcoa_df) <- c("Axis1", "Axis2")
  pcoa_data <- cbind(metadata, pcoa_df)
  
  evals <- pcoa_ord$values$Eigenvalues
  pc1_var <- round(100 * evals[1] / sum(evals), 1)
  pc2_var <- round(100 * evals[2] / sum(evals), 1)
  
  # 4. Perform NMDS Ordination
  capture.output({
    nmds_ord <- try(ordinate(physeq_prop, method = "NMDS", distance = curr_dist, k=2, trymax=50), silent=TRUE)
  })
}

# -------------------------------------------------------------------------
# 3. Differential Abundance with ANCOM-BC2
# -------------------------------------------------------------------------
cat(">>> Section 3: Running ANCOM-BC2 Global Differential Abundance...\n")

library(ANCOMBC)

# Agglomerate to Species
physeq_species <- tax_glom(physeq, taxrank = "Species", NArm = FALSE)

# Clean group names to remove spaces to ensure robust column matching in ANCOM-BC2
sample_data(physeq_species)$Group <- gsub(" ", "", sample_data(physeq_species)$Group)
# Run ANCOM-BC2 (both Global LRT and Pairwise comparisons in one shot)
res_ancom <- ancombc2(
  data = physeq_species,
  fix_formula = "Group",
  group = "Group",
  global = TRUE,
  pairwise = TRUE,
  verbose = FALSE
)

# Extract and format Global results
res_global <- res_ancom$res_global
tax_species <- as.data.frame(tax_table(physeq_species))

diff_df <- res_global %>%
  mutate(
    ASV_ID = taxon,
    Taxon = paste0(tax_species[ASV_ID, "Genus"], "_", tax_species[ASV_ID, "Species"]),
    P_Group = p_val,
    FDR_Group = q_val,
    F_Group = W # W statistic acts as the global test statistic
  ) %>%
  filter(!is.na(P_Group)) %>%
  arrange(P_Group) %>%
  select(ASV_ID, Taxon, F_Group, P_Group, FDR_Group)

write_tsv(diff_df, file.path(STATS_DIR, "differential", "differential_abundance_species.tsv"))
cat("ANCOM-BC2 global differential abundance results saved.\n")

# Taxonomic Composition stacked barplots omitted. Handled dynamically by vis/1.tax.R

cat(">>> Section 5: Pairwise Group-by-Group Differential Abundance...\n")

res_pair <- res_ancom$res_pair

# Save the results
out_file <- file.path(STATS_DIR, "differential", "pairwise_differential_abundance.tsv")
write_tsv(res_pair, out_file)
cat(sprintf("   Pairwise ANCOM-BC2 results saved to: %s\n", out_file))