library(phyloseq)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(htmltools)
library(dplyr)

# 1. Load data
ps <- readRDS("results/03_phyloseq/phyloseq_obj.rds")

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# 2. Function to build and save individual plots
save_taxa_plot <- function(physeq, rank_name, file_out) {
  # Agglomerate and melt
  ps_glom <- tax_glom(physeq, taxrank = rank_name)
  df <- psmelt(ps_glom)
  
  # Standardize Taxon column and handle NAs
  df$Taxon <- as.character(df[[rank_name]])
  df$Taxon[is.na(df$Taxon)] <- "Unclassified"
  
  # Group other taxa outside of top 20 into "Other"
  # Calculate mean relative abundance per taxon to identify top 20
  top_taxa <- df %>%
    group_by(Taxon) %>%
    summarise(MeanAbund = mean(Abundance, na.rm = TRUE)) %>%
    filter(Taxon != "Unclassified") %>%
    arrange(desc(MeanAbund)) %>%
    head(20) %>%
    pull(Taxon)
  
  # Reassign Taxon: if not in top 20 and not Unclassified, rename to "Other"
  df$Taxon <- ifelse(df$Taxon %in% top_taxa | df$Taxon == "Unclassified", df$Taxon, "Other")
  
  # Aggregate abundances by Sample, Taxon, and Sample Variables to avoid multiple segments for "Other"
  group_cols <- c("Sample", "Taxon", sample_variables(physeq))
  group_cols <- intersect(group_cols, colnames(df))
  
  df_grouped <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Define factor levels for Taxon to control legend order (top taxa first, then Unclassified, then Other)
  legend_order <- c(top_taxa)
  if ("Unclassified" %in% df_grouped$Taxon) {
    legend_order <- c(legend_order, "Unclassified")
  }
  if ("Other" %in% df_grouped$Taxon) {
    legend_order <- c(legend_order, "Other")
  }
  
  # Ensure all unique Taxon values are in legend_order
  remaining_taxa <- unique(df_grouped$Taxon)
  remaining_taxa <- remaining_taxa[!remaining_taxa %in% legend_order]
  legend_order <- c(legend_order, remaining_taxa)
  
  df_grouped$Taxon <- factor(df_grouped$Taxon, levels = legend_order)
  
  # Generate a beautiful custom color palette
  num_top <- length(top_taxa)
  if (num_top > 0) {
    top_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(num_top, 12), "Paired"))(num_top)
    names(top_colors) <- top_taxa
  } else {
    top_colors <- c()
  }
  
  custom_colors <- top_colors
  if ("Unclassified" %in% df_grouped$Taxon) {
    custom_colors <- c(custom_colors, "Unclassified" = "#D3D3D3")
  }
  if ("Other" %in% df_grouped$Taxon) {
    custom_colors <- c(custom_colors, "Other" = "#888888")
  }
  
  # Base ggplot
  p <- ggplot(df_grouped, aes(x = Sample, y = Abundance, fill = Taxon, text = paste("Group:", Group))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = custom_colors) +
    theme_minimal() +
    labs(title = paste("Relative Abundance at", rank_name, "Level"), y = "Relative Abundance") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
      legend.position = "right" # Show legend on the right
    )
  
  # Save PDF version
  pdf_out <- sub("\\.html$", ".pdf", file_out)
  pdf_filename <- basename(pdf_out)
  ggsave(pdf_out, plot = p, width = 12, height = 8, device = "pdf", dpi = 300)
  
  # Convert to Plotly and configure high-res download
  p_inter <- ggplotly(p, tooltip = c("fill", "y", "x", "text")) %>%
    layout(
      margin = list(b = 120, r = 350), # 120px bottom margin for rotated sample names, 350px right margin for legend words
      legend = list(orientation = "v", x = 1.02, y = 1)
    ) %>%
    config(
      displaylogo = FALSE, 
      modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "lasso2d", "select2d"),
      toImageButtonOptions = list(
        format = "png", filename = paste0("Taxa_", rank_name), width = 1200, height = 800, scale = 2
      )
    )
  
  # Prepend premium HTML header with 'Download PDF Plot' button
  p_inter_styled <- p_inter %>%
    htmlwidgets::prependContent(
      tags$head(
        tags$style(HTML("
          body {
            font-family: 'Outfit', 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background-color: #fafbfc;
            margin: 0;
            padding: 20px;
            color: #1e293b;
          }
          .report-container {
            max-width: 100%;
            margin: 0 auto;
            background: #ffffff;
            padding: 25px;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.03);
            border: 1px solid #f1f5f9;
            margin-bottom: 20px;
          }
          .header-section {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
            border-bottom: 1px solid #f1f5f9;
            padding-bottom: 15px;
          }
          .title-group h1 {
            font-size: 22px;
            font-weight: 700;
            margin: 0 0 4px 0;
            color: #0f172a;
          }
          .title-group p {
            font-size: 13px;
            color: #64748b;
            margin: 0;
          }
          .btn-download {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: linear-gradient(135deg, #ef4444 0%, #dc2626 100%);
            color: #ffffff !important;
            text-decoration: none !important;
            padding: 8px 16px;
            border-radius: 6px;
            font-weight: 600;
            font-size: 13px;
            box-shadow: 0 4px 10px rgba(239, 68, 68, 0.2);
            transition: all 0.2s ease-in-out;
            border: none;
            cursor: pointer;
          }
          .btn-download:hover {
            transform: translateY(-1px);
            box-shadow: 0 6px 14px rgba(239, 68, 68, 0.3);
            background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%);
          }
          .btn-download:active {
            transform: translateY(1px);
          }
        ")),
        tags$script(HTML("
          document.addEventListener('DOMContentLoaded', function() {
            const resizeObserver = new ResizeObserver(entries => {
              for (let entry of entries) {
                window.dispatchEvent(new Event('resize'));
              }
            });
            resizeObserver.observe(document.body);
          });
        "))
      ),
      tags$div(
        class = "report-container",
        tags$div(
          class = "header-section",
          tags$div(
            class = "title-group",
            tags$h1(paste("Relative Abundance at", rank_name, "Level")),
            tags$p(paste("Interactive visualization of top 20 abundant taxa. All remaining taxa are grouped under 'Other'."))
          ),
          tags$a(
            class = "btn-download",
            href = pdf_filename,
            download = pdf_filename,
            HTML('<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>'),
            "Download PDF Plot"
          )
        )
      )
    )
  
  # Save widget
  saveWidget(p_inter_styled, file_out, selfcontained = TRUE)
}

# 3. Generate robust individual HTML files for ALL ranks
# Add "Kingdom" to this list if you want it, but usually, 16S is all Bacteria.
ranks_to_plot <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

for (rank in ranks_to_plot) {
  tryCatch({
    # Dynamically format the output file name (e.g., 01_taxa_class.html)
    file_out <- paste0("results/interactive_reports/01_taxa_", tolower(rank), ".html")
    save_taxa_plot(ps_rel, rank, file_out)
    message(paste("Successfully generated:", rank))
  }, error = function(e) {
    message(paste("Note: Rank", rank, "not found or failed to plot. Skipping."))
  })
}