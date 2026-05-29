library(ggplot2)
library(plotly)
library(htmlwidgets)
library(htmltools)
library(DT)
library(dplyr)

# ==============================================================================
# PLOT 1: MULTI-GROUP DIFFERENTIAL ABUNDANCE (ANOVA / Kruskal-Wallis)
# ==============================================================================

# 1. Load multi-group data
multi_df <- read.delim("results/04_stats/differential/differential_abundance_genus.tsv", sep = "\t", header = TRUE, check.names = FALSE)

multi_df$Genus <- multi_df$Taxon
multi_df$p_adj <- multi_df$FDR_Group
multi_df$EffectSize <- multi_df$F_Group
p_threshold <- 0.05

multi_df <- multi_df %>%
  filter(!is.na(p_adj) & !is.na(EffectSize)) %>%
  mutate(p_adj = ifelse(p_adj == 0, .Machine$double.xmin, p_adj)) %>%
  mutate(
    Significance = ifelse(p_adj < p_threshold, "Significant Biomarker", "Not Significant"),
    log10_p = -log10(p_adj)
  )

p_multi <- ggplot(multi_df, aes(x = EffectSize, y = log10_p, color = Significance, 
  text = paste("Genus:", Genus, "<br>F-Statistic:", round(EffectSize, 2), "<br>FDR:", signif(p_adj, 3)))) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Significant Biomarker" = "red", "Not Significant" = "grey")) +
  theme_bw() +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "blue", alpha = 0.5) +
  labs(title = "Multi-Group Biomarkers (ANOVA/KW)", x = "F-Statistic (Variance Magnitude)", y = "-Log10(FDR)")

p_multi_inter <- ggplotly(p_multi, tooltip = "text") %>%
  config(displaylogo = FALSE, toImageButtonOptions = list(format="png", filename="MultiGroup_Biomarkers", width=1000, height=700, scale=2))

sig_multi <- multi_df %>% filter(Significance == "Significant Biomarker") %>% arrange(p_adj) %>% select(Genus, EffectSize, P_Group, p_adj)
sig_multi$EffectSize <- round(sig_multi$EffectSize, 2)
sig_multi$P_Group <- signif(sig_multi$P_Group, 3)
sig_multi$p_adj <- signif(sig_multi$p_adj, 3)
colnames(sig_multi) <- c("Genus", "F-Statistic", "P-Value", "FDR")

# Save PDF version
ggsave("results/interactive_reports/04_diff_multi.pdf", plot = p_multi, width = 8, height = 6, device = "pdf", dpi = 300)

multi_view <- tagList(
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
        max-width: 1200px;
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
      .section-block {
        margin-top: 25px;
        border-top: 1px solid #f1f5f9;
        padding-top: 20px;
      }
    "))
  ),
  tags$div(
    class = "report-container",
    tags$div(
      class = "header-section",
      tags$div(
        class = "title-group",
        tags$h3("Global Biomarkers (All Groups)"),
        tags$p("Taxa driving the overall structural variance across all tested groups.")
      ),
      tags$a(
        class = "btn-download",
        href = "04_diff_multi.pdf",
        download = "04_diff_multi.pdf",
        HTML('<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>'),
        "Download PDF Plot"
      )
    ),
    tags$p("Global ANOVA/Kruskal-Wallis analysis to identify significant taxa whose abundances vary across all study groups."),
    p_multi_inter,
    
    tags$div(
      class = "section-block",
      datatable(sig_multi, options = list(pageLength = 5), rownames = FALSE)
    )
  )
)
save_html(multi_view, "results/interactive_reports/04_diff_multi.html")

# ==============================================================================
# PLOT 2: PAIRWISE VOLCANO PLOT (Group vs Group)
# ==============================================================================

# Check if pairwise data exists before trying to plot it
pair_file <- "results/04_stats/differential/pairwise_differential_abundance.tsv"
if (file.exists(pair_file)) {
  pair_df <- read.delim(pair_file, sep = "\t", header = TRUE)
  
  lfc_threshold <- 1.0 # Log2FC > 1 or < -1 (Fold change of 2)
  
  pair_df <- pair_df %>%
    filter(!is.na(p_adj) & !is.na(log2FoldChange)) %>%
    mutate(p_adj = ifelse(p_adj == 0, .Machine$double.xmin, p_adj)) %>%
    mutate(
      Significance = case_when(
        p_adj < p_threshold & log2FoldChange > lfc_threshold ~ "Enriched (Up)",
        p_adj < p_threshold & log2FoldChange < -lfc_threshold ~ "Depleted (Down)",
        TRUE ~ "Not Significant"
      ),
      log10_p = -log10(p_adj)
    )
  
  # Note the 'frame = Comparison' aesthetic. This creates the dropdown animation.
  p_pair <- ggplot(pair_df, aes(x = log2FoldChange, y = log10_p, color = Significance, frame = Comparison,
    text = paste("Genus:", Genus, "<br>Log2FC:", round(log2FoldChange, 2), "<br>FDR:", signif(p_adj, 3)))) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("Enriched (Up)" = "red", "Depleted (Down)" = "blue", "Not Significant" = "grey")) +
    theme_bw() +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
    labs(title = "Pairwise Biomarkers (Volcano Plot)", x = "Log2 Fold Change", y = "-Log10(FDR)")
  
  p_pair_inter <- ggplotly(p_pair, tooltip = "text") %>%
    animation_opts(frame = 0, transition = 0, redraw = TRUE) %>%
    animation_slider(currentvalue = list(prefix = "Comparison: ", font = list(color = "red"))) %>%
    config(displaylogo = FALSE, toImageButtonOptions = list(format="png", filename="Pairwise_Volcano", width=1000, height=700, scale=2))
  
  sig_pair <- pair_df %>% filter(Significance != "Not Significant") %>% arrange(Comparison, p_adj) %>% select(Comparison, Genus, log2FoldChange, p_adj, Significance)
  sig_pair$log2FoldChange <- round(sig_pair$log2FoldChange, 2)
  sig_pair$p_adj <- signif(sig_pair$p_adj, 3)
  
  # Save PDF version (faceted by Comparison so all are visible simultaneously)
  p_pair_static <- p_pair + facet_wrap(~Comparison)
  ggsave("results/interactive_reports/04_diff_pairwise.pdf", plot = p_pair_static, width = 12, height = 8, device = "pdf", dpi = 300)
  
  pair_view <- tagList(
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
          max-width: 1200px;
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
        .section-block {
          margin-top: 25px;
          border-top: 1px solid #f1f5f9;
          padding-top: 20px;
        }
      "))
    ),
    tags$div(
      class = "report-container",
      tags$div(
        class = "header-section",
        tags$div(
          class = "title-group",
          tags$h3("Pairwise Biomarkers (Group-by-Group)"),
          tags$p("Use the slider at the bottom of the plot to toggle between specific group comparisons.")
        ),
        tags$a(
          class = "btn-download",
          href = "04_diff_pairwise.pdf",
          download = "04_diff_pairwise.pdf",
          HTML('<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>'),
          "Download PDF Plot"
        )
      ),
      tags$p("Genera in red are enriched in the first group relative to the second, while blue indicates depletion. The static PDF version is beautifully faceted to display all comparisons simultaneously."),
      p_pair_inter,
      
      tags$div(
        class = "section-block",
        datatable(sig_pair, options = list(pageLength = 5), rownames = FALSE, caption = "Significant Pairwise Shifts")
      )
    )
  )
  save_html(pair_view, "results/interactive_reports/04_diff_pairwise.html")
} else {
  message("Pairwise data not found. Skipping pairwise plot.")
}