#!/usr/bin/env Rscript

# =========================================================================
#          16S rRNA Microbial Co-occurrence Network Analysis
# =========================================================================

cat("Loading required libraries...\n")
library(phyloseq)
library(tidyverse)
library(igraph)
library(ggplot2)
library(cowplot)

# Get parameters from environment variables (with default fallbacks)
PHYLO_DIR <- Sys.getenv("PHYLO_DIR", "./results/03_phyloseq")
NETWORK_DIR <- Sys.getenv("NETWORK_DIR", "./results/05_network")

cat("=========================================================================\n")
cat("Parameters:\n")
cat("Phyloseq Dir: ", PHYLO_DIR, "\n")
cat("Network Dir:  ", NETWORK_DIR, "\n")
cat("=========================================================================\n")

# Create output directory
dir.create(NETWORK_DIR, showWarnings = FALSE, recursive = TRUE)

# Load Phyloseq object
physeq_path <- file.path(PHYLO_DIR, "phyloseq_obj.rds")
if (!file.exists(physeq_path)) {
  stop("ERROR: Phyloseq object not found at ", physeq_path)
}
physeq <- readRDS(physeq_path)

# Extract metadata
metadata <- as(sample_data(physeq), "data.frame")
groups <- unique(metadata$Group)

# Collapse to Species level
cat(">>> Agglomerating to Species level...\n")
physeq_species <- tax_glom(physeq, taxrank = "Species", NArm = FALSE)

cat(">>> Calculating prevalence mask directly from Species object...\n")

# Dynamically handle whether taxa are rows or columns in your matrix
if (taxa_are_rows(physeq_species)) {
  species_counts <- as(otu_table(physeq_species), "matrix")
  prev_counts  <- rowSums(species_counts > 0)
} else {
  species_counts <- as(otu_table(physeq_species), "matrix")
  prev_counts  <- colSums(species_counts > 0)
}

# Define your threshold (e.g., must be present in at least 2 samples)
min_samples_threshold <- 2
prev_mask <- prev_counts >= min_samples_threshold

# CRITICAL FIX: Ensure lengths match identically before pruning
if (length(prev_mask) == ntaxa(physeq_species)) {
  physeq_species_filt <- prune_taxa(prev_mask, physeq_species)
  cat(">>> Successfully filtered out rare species. Retained:", ntaxa(physeq_species_filt), "species.\n")
} else {
  stop("CRITICAL ERROR: Vector lengths still do not match. Check matrix orientations.")
}

# Extract count data and taxonomy mapping
count_matrix <- as(otu_table(physeq_species_filt), "matrix")
if (taxa_are_rows(physeq_species_filt)) {
  count_matrix <- t(count_matrix)
}
tax_species <- as.data.frame(tax_table(physeq_species_filt))

# Dataframes to store network topology comparison metrics
network_metrics <- data.frame(
  Group = character(),
  Nodes = integer(),
  Edges = integer(),
  Density = numeric(),
  AverageDegree = numeric(),
  ClusteringCoefficient = numeric(),
  stringsAsFactors = FALSE
)

# Set up list to store plots
network_plots <- list()

# -------------------------------------------------------------------------
# Build and Compare Networks for Each Group
# -------------------------------------------------------------------------
for (grp in groups) {
  cat("\nProcessing network for Group:", grp, "...\n")
  
  # Subset sample data to current group
  grp_samples <- rownames(metadata[metadata$Group == grp, , drop = FALSE])
  if (length(grp_samples) < 3) {
    cat("Skipping Group", grp, "due to insufficient samples (<3).\n")
    next
  }
  
  # Subset counts
  grp_counts <- count_matrix[grp_samples, , drop = FALSE]
  
  # Calculate Spearman Correlation
  cor_res <- cor(grp_counts, method = "spearman")
  
  # Replace any NA/NaN correlations (due to zero variance/constant values) with 0
  cor_res[is.na(cor_res)] <- 0
  
  # Get p-values for correlations using corr.p approximation
  n <- nrow(grp_counts)
  # Prevent division by zero if cor_res is exactly 1 or -1
  cor_for_t <- cor_res
  cor_for_t[abs(cor_for_t) >= 0.9999] <- 0.9999 * sign(cor_for_t[abs(cor_for_t) >= 0.9999])
  
  t_stat <- cor_for_t * sqrt((n - 2) / (1 - cor_for_t^2))
  p_vals <- 2 * (1 - pt(abs(t_stat), df = n - 2))
  p_vals[is.na(p_vals)] <- 1
  diag(p_vals) <- 1 # diagonal is self-correlation
  
  # Adjust p-values for multiple testing (FDR)
  p_adj <- matrix(p.adjust(as.vector(p_vals), method = "BH"), nrow = nrow(p_vals))
  rownames(p_adj) <- rownames(p_vals)
  colnames(p_adj) <- colnames(p_vals)
  
  # Apply strict filtering: Keep only strong (|r| >= 0.6) and significant (FDR q < 0.05) connections
  cor_adj <- cor_res
  cor_adj[abs(cor_res) < 0.6 | p_adj >= 0.05] <- 0
  diag(cor_adj) <- 0 # remove self loops
  
  # Generate igraph object
  g <- graph_from_adjacency_matrix(cor_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Remove isolated nodes (degree = 0)
  isolated_nodes <- V(g)[degree(g) == 0]
  g <- delete_vertices(g, isolated_nodes)
  
  num_nodes <- vcount(g)
  num_edges <- ecount(g)
  
  cat("Group", grp, "Network contains", num_nodes, "nodes and", num_edges, "edges.\n")
  
  if (num_nodes == 0) {
    cat("No significant co-occurrence relationships found in Group", grp, ". Skipping network plotting.\n")
    next
  }
  
  # Calculate topological properties
  density <- edge_density(g)
  avg_degree <- mean(degree(g))
  clust_coeff <- transitivity(g, type = "global")
  
  # Append to summary metrics table
  network_metrics <- rbind(network_metrics, data.frame(
    Group = grp,
    Nodes = num_nodes,
    Edges = num_edges,
    Density = round(density, 4),
    AverageDegree = round(avg_degree, 2),
    ClusteringCoefficient = round(clust_coeff, 4)
  ))
  
  # Add node labels (Species level) and color by Phylum
  # Retrieve taxonomy for remaining vertices
  node_asvs <- V(g)$name
  node_species <- tax_species[node_asvs, "Species"]
  node_phyla <- tax_species[node_asvs, "Phylum"]
  
  V(g)$label <- node_species
  V(g)$phylum <- node_phyla
  
  # Assign color based on phylum
  phylum_factors <- factor(node_phyla)
  palette <- rainbow(length(levels(phylum_factors)))
  V(g)$color <- palette[as.numeric(phylum_factors)]
  
  # Set edge properties: green for positive correlation, red for negative
  E(g)$color <- ifelse(E(g)$weight > 0, "#2ca02c", "#d62728")
  E(g)$width <- abs(E(g)$weight) * 3
  
  # Save network object
  saveRDS(g, file.path(NETWORK_DIR, paste0("network_", tolower(grp), ".rds")))
  
  # Save individual layout plots as images
  png(file.path(NETWORK_DIR, paste0("network_plot_", tolower(grp), ".png")), width = 800, height = 800)
  set.seed(42) # For reproducible layout
  plot(g, 
       layout = layout_with_fr(g, weights = abs(E(g)$weight)), 
       vertex.size = 6 + degree(g)*0.3,
       vertex.label.color = "black",
       vertex.label.cex = 0.8,
       vertex.label.dist = 1.0,
       main = paste("Microbial Co-occurrence Network:", grp))
  legend("bottomleft", legend = levels(phylum_factors), fill = palette, title = "Phylum", bty = "n")
  dev.off()
}

# -------------------------------------------------------------------------
# Export Network Topology Metrics
# -------------------------------------------------------------------------
write_tsv(network_metrics, file.path(NETWORK_DIR, "network_topology_comparison.tsv"))

# Create a clean visual plot of topological metrics
p_nodes <- ggplot(network_metrics, aes(x = Group, y = Nodes, fill = Group)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_cowplot(12) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Network Nodes", y = "Count") +
  theme(legend.position = "none")

p_edges <- ggplot(network_metrics, aes(x = Group, y = Edges, fill = Group)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_cowplot(12) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Network Edges", y = "Count") +
  theme(legend.position = "none")

p_density <- ggplot(network_metrics, aes(x = Group, y = Density, fill = Group)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_cowplot(12) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Network Density", y = "Density Ratio") +
  theme(legend.position = "none")

p_degree <- ggplot(network_metrics, aes(x = Group, y = AverageDegree, fill = Group)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_cowplot(12) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Average Degree", y = "Mean Connections") +
  theme(legend.position = "none")

combined_metrics_plot <- plot_grid(p_nodes, p_edges, p_density, p_degree, ncol = 2)
ggsave(file.path(NETWORK_DIR, "network_topology_comparison.png"), 
       plot = combined_metrics_plot, width = 10, height = 8, dpi = 300)

cat("\n=========================================================================\n")
cat("Co-occurrence Network Analysis Completed successfully!\n")
cat("Topology Table: ", file.path(NETWORK_DIR, "network_topology_comparison.tsv"), "\n")
cat("Topology Plot:  ", file.path(NETWORK_DIR, "network_topology_comparison.png"), "\n")
cat("=========================================================================\n")
