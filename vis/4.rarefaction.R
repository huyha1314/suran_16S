library(phyloseq)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(htmltools)
library(vegan)
library(dplyr)
library(tidyr)

# 1. Load the pre-compiled phyloseq object
ps <- readRDS("results/03_phyloseq/phyloseq_obj.rds")

meta_df <- as.data.frame(sample_data(ps))
meta_df$SampleID <- rownames(meta_df)

# Create output folder
dir.create("results/interactive_reports", recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# HELPER FUNCTION TO GENERATE RAREFACTION DATA (OPTIMIZED WITH MULTINOMIAL DRAW)
# =========================================================================
get_rarefaction_data <- function(otu_tab, metric = "observed") {
  rarefaction_data <- list()
  
  for (i in 1:nrow(otu_tab)) {
    sample_name <- rownames(otu_tab)[i]
    sample_seq <- otu_tab[i, ]
    sample_seq <- sample_seq[sample_seq > 0] # Remove 0s for speed
    
    n_reads <- sum(sample_seq)
    
    if (n_reads > 0) {
      # Adapt step size to guarantee ~25 smooth steps per sample
      step_size <- max(1000, round(n_reads / 25))
      steps <- unique(c(seq(1, n_reads, by = step_size), n_reads))
      
      if (metric == "observed") {
        # Analytical expected richness using vegan::rarefy (extremely fast)
        vals <- as.numeric(rarefy(sample_seq, steps))
      } else if (metric == "shannon") {
        # Multinomial sampling (blazing fast compared to rrarefy)
        probs <- sample_seq / n_reads
        vals <- sapply(steps, function(d) {
          if (d == n_reads) {
            diversity(sample_seq, index = "shannon")
          } else {
            sub_seq <- rmultinom(1, size = d, prob = probs)[, 1]
            diversity(sub_seq, index = "shannon")
          }
        })
      }
      
      rarefaction_data[[i]] <- data.frame(
        SampleID = sample_name,
        Reads = steps,
        Value = vals
      )
    }
  }
  
  # Bind and join metadata
  rc_df <- bind_rows(rarefaction_data)
  plot_df <- left_join(rc_df, meta_df, by = "SampleID") %>%
    arrange(SampleID, Reads)
  
  return(plot_df)
}

# =========================================================================
# HELPER FUNCTION TO GENERATE INTERACTIVE PLOTLY AND PDF
# =========================================================================
save_rarefaction_report <- function(plot_df, title, y_label, min_depth, file_out, filename) {
  p <- ggplot(plot_df, aes(x = Reads, y = Value, group = SampleID, color = Group, text = SampleID)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    theme_bw() +
    labs(
      title = title,
      x = "Sequencing Depth (Number of Reads)",
      y = y_label
    ) +
    geom_vline(xintercept = min_depth, linetype = "dashed", color = "darkgrey")
  
  # Save PDF version
  pdf_out <- sub("\\.html$", ".pdf", file_out)
  pdf_filename <- basename(pdf_out)
  ggsave(pdf_out, plot = p, width = 10, height = 7, device = "pdf", dpi = 300)
  
  # Interactive Plotly
  p_interactive <- ggplotly(p, tooltip = c("text", "x", "y")) %>%
    layout(
      margin = list(b = 60, r = 180), # 60px bottom margin, 180px right margin for group legend names
      legend = list(orientation = "v", x = 1.02, y = 1)
    ) %>%
    config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "lasso2d", "select2d"),
      toImageButtonOptions = list(
        format = "png",
        filename = filename,
        width = 1000,
        height = 700,
        scale = 2
      )
    )
  
  # Wrap in beautiful premium container with Download button
  p_inter_styled <- p_interactive %>%
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
            tags$h1(title),
            tags$p("Interactive curve showing taxonomic saturation as sequencing depth increases. Dashed line represents minimum sequencing depth.")
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
  
  saveWidget(p_inter_styled, file_out, selfcontained = TRUE)
}

# =========================================================================
# CALCULATE AND BUILD THE PLOTS
# =========================================================================

# 1. ASV Level Raw Counts
otu_asv <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) { otu_asv <- t(otu_asv) }
min_depth_asv <- min(rowSums(otu_asv))

# 2. Genus Level Raw Counts
ps_genus <- tax_glom(ps, taxrank = "Genus")
otu_genus <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) { otu_genus <- t(otu_genus) }
min_depth_genus <- min(rowSums(otu_genus))

cat(">>> Calculating ASV-Level Richness Rarefaction...\n")
df_asv_obs <- get_rarefaction_data(otu_asv, metric = "observed")
save_rarefaction_report(df_asv_obs, "ASV-Level Observed Richness Rarefaction", "Observed ASVs (Richness)", min_depth_asv, "results/interactive_reports/02b_rarefaction_asv_observed.html", "ASV_Observed_Rarefaction")

cat(">>> Calculating Genus-Level Richness Rarefaction...\n")
df_genus_obs <- get_rarefaction_data(otu_genus, metric = "observed")
save_rarefaction_report(df_genus_obs, "Genus-Level Observed Richness Rarefaction", "Observed Genera", min_depth_genus, "results/interactive_reports/02b_rarefaction_genus_observed.html", "Genus_Observed_Rarefaction")

cat(">>> Calculating ASV-Level Shannon Saturation...\n")
df_asv_sha <- get_rarefaction_data(otu_asv, metric = "shannon")
save_rarefaction_report(df_asv_sha, "ASV-Level Shannon Diversity Saturation", "Shannon Diversity Index", min_depth_asv, "results/interactive_reports/02b_rarefaction_asv_shannon.html", "ASV_Shannon_Saturation")

cat(">>> Rarefaction calculations and interactive plotting completed successfully!\n")