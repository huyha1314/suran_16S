library(phyloseq)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(htmltools)
library(DT)
library(dplyr)

# 1. Load the pre-compiled phyloseq object
ps <- readRDS("results/03_phyloseq/phyloseq_obj.rds")

# Transform data to relative abundance for Beta Diversity calculation
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Ensure output directory for plots exists
beta_plot_dir <- "results/04_stats/beta"
dir.create(beta_plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# GENERATE HIGH-RES PNG PLOTS FOR PCoA AND NMDS (ALL METRICS)
# =========================================================================
metrics <- c("bray", "jaccard", "wunifrac", "unifrac")
metrics_file_names <- c("bray_curtis", "jaccard", "weighted_unifrac", "unweighted_unifrac")
metric_names <- c("Bray-Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac")

has_tree <- !is.null(phy_tree(ps_rel, errorIfNULL = FALSE))

for (i in 1:length(metrics)) {
  m <- metrics[i]
  m_file <- metrics_file_names[i]
  m_name <- metric_names[i]
  
  # Skip UniFrac metrics if no phylogenetic tree is available
  if (m %in% c("wunifrac", "unifrac") && !has_tree) {
    cat(sprintf("   Skipping %s: no phylogenetic tree available (long-read mode).\n", m_name))
    next
  }
  
  cat(sprintf(">>> Generating high-res Beta Diversity plots for %s...\n", m_name))
  
  # 1. PCoA Ordination
  ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = m)
  p_pcoa <- plot_ordination(ps_rel, ord_pcoa, color = "Group") +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Group), type = "norm", linetype = 2, alpha = 0.4) +
    theme_bw() +
    labs(
      title = sprintf("Beta Diversity PCoA (%s)", m_name),
      color = "Group"
    )
  ggsave(file.path(beta_plot_dir, sprintf("pcoa_%s.png", m_file)), plot = p_pcoa, width = 8, height = 6, dpi = 300)
  
  # 2. NMDS Ordination
  capture.output({
    ord_nmds <- try(ordinate(ps_rel, method = "NMDS", distance = m, k = 2, trymax = 50), silent = TRUE)
  })
  
  if (!inherits(ord_nmds, "try-error")) {
    p_nmds <- plot_ordination(ps_rel, ord_nmds, color = "Group") +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(aes(group = Group), type = "norm", linetype = 2, alpha = 0.4) +
      theme_bw() +
      labs(
        title = sprintf("Beta Diversity NMDS (%s)", m_name),
        subtitle = sprintf("Stress = %s", round(ord_nmds$stress, 4)),
        color = "Group"
      )
    ggsave(file.path(beta_plot_dir, sprintf("nmds_%s.png", m_file)), plot = p_nmds, width = 8, height = 6, dpi = 300)
  } else {
    cat(sprintf("   Warning: NMDS ordination failed to converge for %s.\n", m_name))
  }
}

# =========================================================================
# BUILD THE MAIN BRAY-CURTIS INTERACTIVE PLOT FOR THE HTML REPORT
# =========================================================================
cat(">>> Building primary Bray-Curtis interactive ordination...\n")
ord_bray <- ordinate(ps_rel, method = "PCoA", distance = "bray")

p_beta <- plot_ordination(ps_rel, ord_bray, color = "Group") +
  geom_point(aes(text = paste("Sample:", rownames(sample_data(ps_rel)))), size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Group), type = "norm", linetype = 2, alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Beta Diversity: PCoA (Bray-Curtis)", 
    color = "Group"
  )

p_beta_interactive <- ggplotly(p_beta, tooltip = c("text", "color")) %>%
  config(
    displaylogo = FALSE, 
    modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "lasso2d", "select2d"),
    toImageButtonOptions = list(
      format = "png",
      filename = "Beta_Diversity_PCoA_Bray_Curtis",
      width = 1000,
      height = 700,
      scale = 2
    )
  )

# =========================================================================
# COMPILING GLOBAL STATISTICAL TABLE
# =========================================================================
summary_list <- list()

for (i in 1:length(metrics_file_names)) {
  m_file <- metrics_file_names[i]
  m_name <- metric_names[i]
  
  perm_file <- file.path("results/04_stats/beta", paste0("permanova_", m_file, ".tsv"))
  disp_file <- file.path("results/04_stats/beta", paste0("permdisp_", m_file, ".tsv"))
  
  perm_p <- NA; perm_r2 <- NA; disp_p <- NA
  
  if (file.exists(perm_file)) {
    perm_df <- read.delim(perm_file, sep="\t", row.names=1, check.names=FALSE)
    if (nrow(perm_df) > 0) {
      perm_p <- perm_df[1, "Pr(>F)"]
      perm_r2 <- perm_df[1, "R2"]
    }
  }
  
  if (file.exists(disp_file)) {
    disp_df <- read.delim(disp_file, sep="\t", row.names=1, check.names=FALSE)
    if (nrow(disp_df) > 0) {
      disp_p <- disp_df[1, "Pr(>F)"]
    }
  }
  
  if (!is.na(perm_p)) {
    summary_list[[i]] <- data.frame(
      DistanceMetric = m_name,
      PERMANOVA_R2 = round(perm_r2, 3),
      PERMANOVA_P = perm_p,
      PERMDISP_P = disp_p,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- bind_rows(summary_list)

if (nrow(summary_df) > 0) {
  summary_df <- summary_df %>%
    mutate(
      PERMANOVA_Sig = ifelse(PERMANOVA_P < 0.05, "Yes (*)", "No"),
      PERMDISP_Sig = ifelse(is.na(PERMDISP_P), "Pending", 
                            ifelse(PERMDISP_P < 0.05, "Unequal Variance (Warning)", "Equal Variance (Good)"))
    )
  
  colnames(summary_df) <- c("Distance Metric", "PERMANOVA R-Squared", "PERMANOVA P-Value", "PERMDISP P-Value", "PERMANOVA Significant", "PERMDISP Assumption")
} else {
  summary_df <- data.frame(Message = "Data is still generating...")
}

summary_table <- datatable(
  summary_df,
  options = list(dom = 't', ordering = FALSE, pageLength = 10),
  rownames = FALSE,
  caption = "Table 1: Beta Diversity Summary across all Distance Metrics"
)

# =========================================================================
# COMBINE EVERYTHING INTO HTML WIDGET
# =========================================================================
# Save PDF version
pdf_bray_out <- "results/interactive_reports/03_beta_diversity_pcoa.pdf"
ggsave(pdf_bray_out, plot = p_beta, width = 8, height = 6, device = "pdf", dpi = 300)

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
        tags$h3("Beta Diversity (Between-Sample Dissimilarity)"),
        tags$p("Principal Coordinates Analysis (PCoA) illustrating the structural similarity of the microbial communities.")
      ),
      tags$a(
        class = "btn-download",
        href = "03_beta_diversity_pcoa.pdf",
        download = "03_beta_diversity_pcoa.pdf",
        HTML('<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>'),
        "Download PDF Plot"
      )
    ),
    tags$p("The Principal Coordinates Analysis (PCoA) plot illustrates the structural similarity of the microbial communities. Samples that are closer together have a more similar microbiome composition. The dotted lines represent 95% confidence intervals for each group."),
    p_beta_interactive,
    
    tags$div(
      class = "section-block",
      tags$h4("Global Statistical Summary"),
      tags$p("The table below summarizes the PERMANOVA and PERMDISP results across all calculated distance metrics. R-Squared indicates how much of the variation is explained by the grouping. A significant PERMANOVA indicates structural differences, while a non-significant PERMDISP validates that the result is not just an artifact of unequal variance."),
      summary_table
    )
  )
)

save_html(combined_view, "results/interactive_reports/03_beta_diversity_with_stats.html")
cat(">>> Beta diversity visualization with high-res PNG and PDF plots completed successfully!\n")