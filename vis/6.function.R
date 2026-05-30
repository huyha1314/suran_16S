library(ggplot2)
library(plotly)
library(htmlwidgets)
library(htmltools)
library(DT)
library(readr)
library(dplyr)
library(tidyr)
library(matrixStats) # For rowVars

# 1. Load the unstratified pathway abundance from PICRUSt2
# This file contains the raw abundance of functional pathways predicted for each sample
picrust_path <- "results/06_picrust2/pathways_out/path_abun_unstrat.tsv"
if (!file.exists(picrust_path)) {
  message("PICRUSt2 output not found (expected in longread mode). Skipping functional visualization.")
  quit(save = "no", status = 0)
}
pathways_raw <- read_tsv(picrust_path, show_col_types = FALSE)

# 2. Filter for the Top Most Variable Pathways
# We want to remove static "housekeeping" pathways and highlight the ones driving functional changes.
pathway_matrix <- as.matrix(pathways_raw[, -1]) # Remove the pathway name column for math
rownames(pathway_matrix) <- pathways_raw$pathway

# Calculate the variance for each pathway across all samples
pathway_variances <- rowVars(pathway_matrix)

# Select the top 40 most variable pathways (you can adjust this number)
top_n <- 40
top_indices <- order(pathway_variances, decreasing = TRUE)[1:top_n]
top_pathways_matrix <- pathway_matrix[top_indices, ]

# 3. Data Transformation for the Heatmap
# We convert the matrix back to a long-format dataframe for ggplot
top_pathways_df <- as.data.frame(top_pathways_matrix)
top_pathways_df$Pathway <- rownames(top_pathways_df)

path_long <- top_pathways_df %>%
  pivot_longer(cols = -Pathway, names_to = "SampleID", values_to = "Abundance")

# We use log10(x + 1) transformation to normalize massive functional counts for visual clarity
path_long$Log10_Abundance <- log10(path_long$Abundance + 1)

# Ensure Pathway names are factors ordered by clustering (optional, but makes heatmap cleaner)
# Here we order them alphabetically for simplicity, but you could use hclust here.
path_long$Pathway <- factor(path_long$Pathway, levels = rev(sort(unique(path_long$Pathway))))

# 4. Create the base ggplot Heatmap
p_heat <- ggplot(path_long, aes(
  x = SampleID, 
  y = Pathway, 
  fill = Log10_Abundance, 
  text = paste(
    "Sample:", SampleID,
    "<br>Pathway:", Pathway,
    "<br>Raw Abundance:", round(Abundance, 0)
  )
)) +
  geom_tile(color = "white") +
  # Using the Viridis color scale which is colorblind friendly and standard for heatmaps
  scale_fill_viridis_c(name = "Log10\nAbundance") + 
  theme_minimal() +
  labs(title = "Predicted Functional Pathways (Top 40 Variable)",
       x = "Samples",
       y = "MetaCyc Pathways") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 8)
  )

# 5. Convert to interactive Plotly
p_heat_interactive <- ggplotly(p_heat, tooltip = "text") %>%
  config(
    displaylogo = FALSE, 
    modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "lasso2d", "select2d"),
    toImageButtonOptions = list(
      format = "png",
      filename = "Functional_Pathways_Heatmap",
      width = 1200,
      height = 900,
      scale = 2
    )
  )

# 6. Save PDF version of heatmap
ggsave("results/interactive_reports/05_functional_pathways.pdf", plot = p_heat, width = 12, height = 10, device = "pdf", dpi = 300)

# 7. Combine Plot and Text into one HTML widget
combined_view <- tagList(
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
        max-width: 1300px;
        margin: 0 auto;
        background: #ffffff;
        padding: 25px;
        border-radius: 12px;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.03);
        border: 1px solid #f1f5f9;
        margin-bottom: 25px;
      }
      .header-section {
        display: flex;
        justify-content: space-between;
        align-items: center;
        margin-bottom: 20px;
        border-bottom: 1px solid #f1f5f9;
        padding-bottom: 15px;
      }
      .title-group h3 {
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
    "))
  ),
  tags$div(
    class = "report-container",
    tags$div(
      class = "header-section",
      tags$div(
        class = "title-group",
        tags$h3("Functional Potential (PICRUSt2)"),
        tags$p("Predicted metabolic pathways present in the microbiome based on 16S taxonomic profiles.")
      ),
      tags$a(
        class = "btn-download",
        href = "05_functional_pathways.pdf",
        download = "05_functional_pathways.pdf",
        HTML('<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>'),
        "Download PDF Plot"
      )
    ),
    tags$p("Based on the 16S taxonomic profiles, PICRUSt2 predicts the metabolic pathways present in the microbiome. The heatmap below displays the top 40 most variable MetaCyc pathways across all samples. Darker colors (purple) indicate low abundance, while brighter colors (yellow) indicate high abundance. Hover over a cell to see the specific raw abundance value."),
    p_heat_interactive
  )
)

# 8. Save as a self-contained HTML widget
save_html(combined_view, "results/interactive_reports/05_functional_pathways.html")