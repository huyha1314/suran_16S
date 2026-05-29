library(ggplot2)
library(plotly)
library(htmlwidgets)
library(DT)         # Adds interactive tables
library(htmltools)  # Allows combining multiple HTML elements
library(dplyr)      # Added for statistical data manipulation

# 1. Load the pre-calculated alpha diversity table
alpha_df <- read.delim("results/04_stats/alpha/alpha_diversity_table.tsv", sep = "\t", header = TRUE)

# Ensure Group is a character/factor so we can loop through it
alpha_df$Group <- as.character(alpha_df$Group)

# 2. Create the base ggplot - define tooltip aesthetic only on geom_jitter to avoid squishing geom_boxplot
p_alpha <- ggplot(alpha_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) + 
  geom_jitter(aes(text = paste0("Sample: ", SampleID, "<br>Group: ", Group, "<br>Shannon: ", round(Shannon, 3))), width = 0.2, size = 1.5, color = "black", alpha = 0.6) +
  theme_bw() +
  labs(title = "Alpha Diversity (Shannon Index)", y = "Shannon Diversity") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold")
  )

# Save PDF version
pdf_out <- "results/interactive_reports/02_alpha_diversity_with_table.pdf"
pdf_filename <- basename(pdf_out)
ggsave(pdf_out, plot = p_alpha, width = 8, height = 6, device = "pdf", dpi = 300)

# 3. Convert to interactive Plotly
p_alpha_interactive <- ggplotly(p_alpha, tooltip = "text") %>%
  layout(
    margin = list(b = 60, l = 60, r = 30, t = 50)
  ) %>%
  config(
    displaylogo = FALSE, 
    modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "lasso2d", "select2d"),
    toImageButtonOptions = list(
      format = "png",
      filename = "Alpha_Diversity_Shannon",
      width = 1000,
      height = 700,
      scale = 2
    )
  )

# 4. Calculate Pairwise Statistical Comparisons (Wilcoxon Rank-Sum)
# Get all unique combinations of groups
groups <- unique(alpha_df$Group)
group_pairs <- combn(groups, 2, simplify = FALSE)

# Loop through each pair and calculate the p-value
stats_list <- lapply(group_pairs, function(pair) {
  g1 <- pair[1]
  g2 <- pair[2]
  
  # Subset data for just these two groups
  data_sub <- alpha_df[alpha_df$Group %in% c(g1, g2), ]
  
  # Run the test safely
  res <- tryCatch({
    wilcox.test(Shannon ~ Group, data = data_sub, exact = FALSE)$p.value
  }, error = function(e) NA)
  
  data.frame(
    Comparison = paste0(g1, " vs ", g2),
    P_Value = res
  )
})

# Combine results and apply FDR correction
stats_df <- bind_rows(stats_list) %>%
  filter(!is.na(P_Value)) %>%
  mutate(
    Adj_P_Value = p.adjust(P_Value, method = "BH"), # Benjamini-Hochberg correction
    Significance = case_when(
      Adj_P_Value < 0.001 ~ "***",
      Adj_P_Value < 0.01  ~ "**",
      Adj_P_Value < 0.05  ~ "*",
      TRUE                ~ "ns" # Not significant
    )
  )

# Save pairwise Wilcoxon stats as a TSV file
write.table(stats_df, "results/04_stats/alpha/pairwise_wilcoxon_shannon.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Format numbers for a clean table
stats_df$P_Value <- signif(stats_df$P_Value, 3)
stats_df$Adj_P_Value <- signif(stats_df$Adj_P_Value, 3)

# Create the DT Stats Table
stats_table <- datatable(
  stats_df,
  options = list(dom = 't', ordering = FALSE), # 't' hides search bar since it's a small stats table
  rownames = FALSE,
  caption = "Table 1: Pairwise Statistical Significance (Wilcoxon Rank-Sum with FDR Correction)"
)

# 5. Create an interactive Data Table for the raw customer data
alpha_df_clean <- alpha_df
alpha_df_clean$Shannon <- round(alpha_df_clean$Shannon, 3)
alpha_df_clean$Chao1 <- round(alpha_df_clean$Chao1, 3)
if("se.chao1" %in% colnames(alpha_df_clean)) {
  alpha_df_clean$se.chao1 <- round(alpha_df_clean$se.chao1, 3)
}

raw_table <- datatable(
  alpha_df_clean,
  options = list(
    pageLength = 5,         
    scrollX = TRUE,         
    dom = 'Bfrtip'          
  ),
  rownames = FALSE,
  caption = "Table 2: Raw Alpha Diversity Metrics per Sample"
)

# 6. Prepend premium HTML header with 'Download PDF Plot' button
styled_header <- tags$span(
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
        tags$h1("Alpha Diversity (Within-Sample Diversity)"),
        tags$p("The boxplot below shows the distribution of the Shannon Index across groups. Hover over points to identify specific samples.")
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

combined_view <- tagList(
  styled_header,
  p_alpha_interactive,
  tags$br(),
  tags$hr(),
  tags$h4("Statistical Comparisons"),
  tags$p("The table below details the statistical difference between groups. An adjusted p-value < 0.05 is considered significant (*)."),
  stats_table,
  tags$br(),
  tags$hr(),
  tags$h4("Raw Metric Data"),
  raw_table
)

# 7. Save as a single, self-contained HTML file
save_html(combined_view, "results/interactive_reports/02_alpha_diversity_with_table.html")
cat(">>> Alpha diversity visualization with high-res PDF and download button completed successfully!\n")