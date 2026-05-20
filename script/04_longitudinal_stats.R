#!/usr/bin/env Rscript

# =========================================================================
#       16S rRNA Longitudinal Statistical Analysis & Visualization
# =========================================================================

cat("Loading required libraries...\n")
library(phyloseq)
library(vegan)
library(tidyverse)
library(lme4)
library(lmerTest)
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

# Validate metadata columns for longitudinal LMM
required_cols <- c("Group", "Day", "SubjectID")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("ERROR: Metadata is missing required columns for longitudinal analysis: ", 
       paste(missing_cols, collapse = ", "), 
       "\nEnsure 'Group', 'Day', and 'SubjectID' are populated.")
}

# Ensure factors are set correctly
metadata$Group <- factor(metadata$Group)
metadata$Day <- factor(metadata$Day)
metadata$SubjectID <- factor(metadata$SubjectID)

# -------------------------------------------------------------------------
# 1. Longitudinal Alpha Diversity Analysis (LMM)
# -------------------------------------------------------------------------
cat(">>> Section 1: Running Alpha Diversity LMM Analysis...\n")

# Compute alpha diversity metrics
alpha_meas <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Chao1"))
alpha_data <- cbind(metadata, alpha_meas)
write.table(alpha_data, file.path(STATS_DIR, "alpha", "alpha_diversity_table.tsv"), 
            sep = "\t", quote = FALSE, col.names = NA)

# Loop through Shannon and Observed metrics
metrics <- c("Shannon", "Observed")
for (metric in metrics) {
  cat("Modeling longitudinal", metric, "with LMM...\n")
  
  # Fit Linear Mixed-Effects Model
  # Fixed effects: Group, Day, and Group:Day interaction
  # Random effect: (1 | SubjectID) - accounts for repeated measures within subjects
  formula_str <- paste(metric, "~ Group * Day + (1 | SubjectID)")
  lmm_model <- lmer(as.formula(formula_str), data = alpha_data)
  
  # Extract ANOVA table (Type III with Satterthwaite approximation)
  anova_res <- anova(lmm_model)
  write.table(as.data.frame(anova_res), 
              file.path(STATS_DIR, "alpha", paste0("lmm_anova_", tolower(metric), ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # Extract summary (specific coefficients)
  summary_res <- summary(lmm_model)
  write.table(as.data.frame(summary_res$coefficients), 
              file.path(STATS_DIR, "alpha", paste0("lmm_coefficients_", tolower(metric), ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # Visualizations
  # Plot 1: Boxplot across groups and days
  p1 <- ggplot(alpha_data, aes_string(x = "Day", y = metric, fill = "Group")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), size = 1.5, alpha = 0.8) +
    theme_cowplot(12) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste(metric, "by Group & Day"), x = "Collection Day", y = metric)
  
  # Plot 2: Volatility / Spaghetti Plot (tracking subjects over time)
  p2 <- ggplot(alpha_data, aes_string(x = "Day", y = metric, group = "SubjectID", color = "Group")) +
    geom_line(alpha = 0.5, size = 1) +
    geom_point(size = 2) +
    theme_cowplot(12) +
    scale_color_brewer(palette = "Set2") +
    facet_wrap(~Group) +
    labs(title = paste(metric, "Subject Volatility"), x = "Collection Day", y = metric) +
    theme(legend.position = "none")
  
  # Combine plots
  combined_p <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1))
  ggsave(file.path(STATS_DIR, "alpha", paste0("alpha_plot_", tolower(metric), ".png")), 
         plot = combined_p, width = 10, height = 8, dpi = 300)
}

# -------------------------------------------------------------------------
# 1.5 Rarefaction Curve Analysis
# -------------------------------------------------------------------------
cat(">>> Section 1.5: Generating Rarefaction Curves...\n")

# Extract raw counts (non-normalized)
counts_raw <- as.matrix(otu_table(physeq))
if (taxa_are_rows(physeq)) {
  counts_raw <- t(counts_raw)
}

# Find sample depths and step sizes
depths <- rowSums(counts_raw)
max_depth <- max(depths)
step_size <- max(100, floor(max_depth / 50))

# Subsample at multiple intervals and calculate observed species
rarefaction_data <- list()
for (i in 1:nrow(counts_raw)) {
  sample_name <- rownames(counts_raw)[i]
  sample_grp <- metadata[sample_name, "Group"]
  sample_day <- metadata[sample_name, "Day"]
  sample_seq <- counts_raw[i, ]
  sample_seq <- sample_seq[sample_seq > 0]
  
  n_reads <- sum(sample_seq)
  steps <- unique(c(seq(1, n_reads, by = step_size), n_reads))
  
  # Calculate expected richness at each step using vegan::rarefy
  richness <- rarefy(sample_seq, steps)
  
  rarefaction_data[[i]] <- data.frame(
    SampleID = sample_name,
    Group = sample_grp,
    Day = sample_day,
    SubsampleDepth = steps,
    ObservedRichness = richness
  )
}
rarefaction_df <- do.call(rbind, rarefaction_data)

# Plot Rarefaction Curve
p_rare <- ggplot(rarefaction_df, aes(x = SubsampleDepth, y = ObservedRichness, group = SampleID, color = Group)) +
  geom_line(alpha = 0.6, size = 0.8) +
  theme_cowplot(12) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Rarefaction Curves by Group",
       x = "Sequencing Depth (Reads)",
       y = "Observed ASVs (Richness)")

ggsave(file.path(STATS_DIR, "alpha", "rarefaction_curve.png"), plot = p_rare, width = 8, height = 6, dpi = 300)


# -------------------------------------------------------------------------
# 2. Longitudinal Beta Diversity (Stratified PERMANOVA & PCoA)
# -------------------------------------------------------------------------
cat(">>> Section 2: Running Beta Diversity Stratified PERMANOVA...\n")

# Normalize counts using relative abundance (Proportions)
physeq_prop <- transform_sample_counts(physeq, function(x) x / sum(x))

# Calculate Bray-Curtis and Jaccard distances
bray_dist <- phyloseq::distance(physeq_prop, method = "bray")
jaccard_dist <- phyloseq::distance(physeq_prop, method = "jaccard")

dists <- list("Bray-Curtis" = bray_dist, "Jaccard" = jaccard_dist)

for (dist_name in names(dists)) {
  curr_dist <- dists[[dist_name]]
  clean_name <- tolower(gsub("-", "_", dist_name))
  
  cat("Analyzing", dist_name, "beta diversity...\n")
  
  # Run Stratified PERMANOVA using adonis2
  # Critical Innovation: to account for repeated measures, we restrict permutations 
  # WITHIN subjects by setting strata = SubjectID.
  # Otherwise, temporal correlation within individuals violates independent permutation assumptions.
  permanova_res <- adonis2(curr_dist ~ Group * Day, 
                           data = metadata, 
                           strata = metadata$SubjectID,
                           permutations = 999)
  
  write.table(as.data.frame(permanova_res), 
              file.path(STATS_DIR, "beta", paste0("permanova_stratified_", clean_name, ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # Perform PCoA Ordination
  pcoa_ord <- ordinate(physeq_prop, method = "MDS", distance = curr_dist)
  
  # Build trajectory visualization (connecting the same subject across days)
  pcoa_df <- data.frame(pcoa_ord$vectors[, 1:2])
  colnames(pcoa_df) <- c("Axis1", "Axis2")
  pcoa_data <- cbind(metadata, pcoa_df)
  
  # Calculate eigenvalues for axes percentage
  evals <- pcoa_ord$values$Eigenvalues
  pc1_var <- round(100 * evals[1] / sum(evals), 1)
  pc2_var <- round(100 * evals[2] / sum(evals), 1)
  
  # Plot PCoA with trajectories
  p_beta <- ggplot(pcoa_data, aes(x = Axis1, y = Axis2, color = Group)) +
    # Draw paths connecting subjects sequentially across Days
    geom_path(aes(group = SubjectID), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), alpha = 0.4, size = 0.8) +
    geom_point(aes(shape = Day), size = 3.5, alpha = 0.9) +
    theme_cowplot(12) +
    scale_color_brewer(palette = "Set2") +
    labs(title = paste0("PCoA based on ", dist_name),
         x = paste0("PC1 (", pc1_var, "%)"),
         y = paste0("PC2 (", pc2_var, "%)")) +
    theme(legend.box = "horizontal")
  
  ggsave(file.path(STATS_DIR, "beta", paste0("pcoa_", clean_name, ".png")), 
         plot = p_beta, width = 8, height = 6, dpi = 300)
}

# -------------------------------------------------------------------------
# 3. Compositional Differential Abundance with CLR-LMM
# -------------------------------------------------------------------------
cat(">>> Section 3: Running Compositional LMM Differential Abundance...\n")

# Collapse to Genus level
physeq_genus <- tax_glom(physeq, taxrank = "Genus", NArm = FALSE)

# Filter low prevalence genera to increase power (keep if present in at least 20% of samples)
prev_threshold <- 0.20
prev_mask <- apply(otu_table(physeq_genus) > 0, 2, sum) >= (prev_threshold * nsamples(physeq_genus))
physeq_genus_filt <- prune_taxa(prev_mask, physeq_genus)

cat("Genera kept for testing after prevalence filter (>=20%):", ntaxa(physeq_genus_filt), "\n")

# Extract count matrix and add pseudo-count for CLR
count_matrix <- as.matrix(otu_table(physeq_genus_filt))
if (taxa_are_rows(physeq_genus_filt)) {
  count_matrix <- t(count_matrix)
}
count_matrix <- count_matrix + 1  # Add pseudocount of 1 to handle zero counts in log

# Perform CLR transformation (Centered Log-Ratio)
# CLR(x) = log(x) - mean(log(x))
clr_matrix <- t(apply(count_matrix, 1, function(x) {
  log_x <- log(x)
  return(log_x - mean(log_x))
}))

# Extract taxonomy for annotation
tax_genus <- as.data.frame(tax_table(physeq_genus_filt))

# List to store results
lmm_diff_results <- list()

for (i in 1:ncol(clr_matrix)) {
  asv_name <- colnames(clr_matrix)[i]
  genus_name <- tax_genus[asv_name, "Genus"]
  family_name <- tax_genus[asv_name, "Family"]
  full_label <- paste0(family_name, "_", genus_name)
  
  # Prepare model data
  model_df <- data.frame(
    CLR_abundance = clr_matrix[, i],
    Group = metadata$Group,
    Day = metadata$Day,
    SubjectID = metadata$SubjectID
  )
  
  # Fit LMM: test Group, Day, and Group:Day interaction
  # We wrap this in tryCatch in case a rare taxon model fails to converge
  tryCatch({
    fit <- lmer(CLR_abundance ~ Group * Day + (1 | SubjectID), data = model_df)
    an_res <- anova(fit)
    
    # Extract p-values for Group, Day, and Group:Day
    p_group <- an_res["Group", "Pr(>F)"]
    p_day <- an_res["Day", "Pr(>F)"]
    p_inter <- an_res["Group:Day", "Pr(>F)"]
    
    # Extract F-statistics
    f_group <- an_res["Group", "F value"]
    f_day <- an_res["Day", "F value"]
    f_inter <- an_res["Group:Day", "F value"]
    
    lmm_diff_results[[i]] <- data.frame(
      ASV_ID = asv_name,
      Taxon = full_label,
      F_Group = f_group,
      P_Group = p_group,
      F_Day = f_day,
      P_Day = p_day,
      F_Interaction = f_inter,
      P_Interaction = p_inter
    )
  }, error = function(e) {
    # If convergence fails, skip
  })
}

# Bind and adjust for multiple testing (FDR correction)
if (length(lmm_diff_results) > 0) {
  diff_df <- do.call(rbind, lmm_diff_results)
  
  # Benjamini-Hochberg FDR correction
  diff_df$FDR_Group <- p.adjust(diff_df$P_Group, method = "BH")
  diff_df$FDR_Day <- p.adjust(diff_df$P_Day, method = "BH")
  diff_df$FDR_Interaction <- p.adjust(diff_df$P_Interaction, method = "BH")
  
  # Sort by Group significance
  diff_df <- diff_df %>% arrange(P_Group)
  
  write_tsv(diff_df, file.path(STATS_DIR, "differential", "lmm_differential_abundance_genus.tsv"))
  cat("Differential abundance results saved.\n")
  
  # Visualizing the top significant Genus
  top_genus <- diff_df$ASV_ID[1]
  top_label <- diff_df$Taxon[1]
  
  if (!is.na(top_genus)) {
    plot_df <- data.frame(
      CLR_abundance = clr_matrix[, top_genus],
      Group = metadata$Group,
      Day = metadata$Day,
      SubjectID = metadata$SubjectID
    )
    
    p_diff <- ggplot(plot_df, aes(x = Day, y = CLR_abundance, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.7) +
      theme_cowplot(12) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste("Abundance of Top Taxon:", top_label),
           x = "Collection Day",
           y = "CLR-Transformed Abundance")
    
    ggsave(file.path(STATS_DIR, "differential", "top_significant_genus_plot.png"), 
           plot = p_diff, width = 8, height = 6, dpi = 300)
  }
} else {
  cat("WARNING: No LMM models converged successfully.\n")
}

# -------------------------------------------------------------------------
# 4. Taxonomic Composition Volatility
# -------------------------------------------------------------------------
cat(">>> Section 4: Generating Taxonomic Volatility Charts...\n")

# Agglomerate to Phylum level
physeq_phylum <- tax_glom(physeq, taxrank = "Phylum")
phylum_df <- psmelt(transform_sample_counts(physeq_phylum, function(x) x / sum(x)))

# Filter to top 5 Phyla, group others
top5_phyla <- phylum_df %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(Abundance)) %>%
  arrange(desc(MeanAbund)) %>%
  slice(1:5) %>%
  pull(Phylum)

phylum_df <- phylum_df %>%
  mutate(PhylumGroup = ifelse(Phylum %in% top5_phyla, Phylum, "Other Phyla"))

# Volatility Plot for Phyla
p_phylum_time <- ggplot(phylum_df, aes(x = Day, y = Abundance, fill = PhylumGroup)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.6) +
  theme_cowplot(12) +
  scale_fill_brewer(palette = "Spectral") +
  facet_wrap(~Group) +
  labs(title = "Phylum Relative Abundance Volatility", x = "Collection Day", y = "Mean Relative Abundance", fill = "Phylum")

ggsave(file.path(STATS_DIR, "taxa_volatility_phylum.png"), 
       plot = p_phylum_time, width = 10, height = 6, dpi = 300)

# Genus relative abundance stacked barplot
cat(">>> Section 4.5: Generating Genus-level Stacked Barplots...\n")

# Use previously collapsed genus phyloseq object, calculate relative abundances
physeq_genus_rel <- transform_sample_counts(physeq_genus, function(x) x / sum(x))
genus_df <- psmelt(physeq_genus_rel)

# Filter to top 10 Genera, group others as "Other Genera"
top10_genera <- genus_df %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(Abundance)) %>%
  arrange(desc(MeanAbund)) %>%
  slice(1:10) %>%
  pull(Genus)

genus_df <- genus_df %>%
  mutate(GenusGroup = ifelse(Genus %in% top10_genera, Genus, "Other Genera"))

# Ensure GenusGroup displays cleanly
genus_df$GenusGroup <- factor(genus_df$GenusGroup, levels = c(top10_genera, "Other Genera"))

# Plot Genus Stacked Barplot faceted by Group
p_genus_bar <- ggplot(genus_df, aes(x = Day, y = Abundance, fill = GenusGroup)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.6) +
  theme_cowplot(12) +
  scale_fill_brewer(palette = "Paired") +
  facet_wrap(~Group) +
  labs(title = "Genus Relative Abundance Volatility", 
       x = "Collection Day", 
       y = "Mean Relative Abundance", 
       fill = "Genus")

ggsave(file.path(STATS_DIR, "taxa_barplot_genus.png"), 
       plot = p_genus_bar, width = 12, height = 7, dpi = 300)


cat("=========================================================================\n")
cat("Longitudinal Statistical Analysis Completed Successfully!\n")
cat("All plots and tables exported to: ", STATS_DIR, "\n")
cat("=========================================================================\n")
