#!/usr/bin/env Rscript
# =============================================================================
# TCR Diversity Analysis
#
# Analyzes T cell receptor diversity metrics across sample types (PBMC, Lung, Tumor):
#   - Shannon Entropy
#   - Normalized Shannon Entropy
#   - Clonality Index
#   - Relative Clonotype Richness
#   - Total productive TCRs
#   - Unique TCRs
#   - Jaccard Index (sample similarity)
#   - Horn Index (sample similarity)
#
# Prerequisite: perl all.sample.TCR.clone.matrix.pl
# Input:  data/TCR/TCR_Statistics.all.txt
#         data/TCR/compare_TCR_Statistics.all.txt
# Output: results/TCR_analysis/allsingleT.TCR_BNT.Statistics.pdf
# =============================================================================

source("scripts/00_load_packages.R")

# Read TCR diversity statistics
data <- read.table("data/TCR/TCR_Statistics.all.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
data$type <- factor(data$type, levels = c("PBMC", "Lung", "Tumor"))

comparisons <- list(c("PBMC", "Lung"), c("Lung", "Tumor"), c("PBMC", "Tumor"))
pal <- c("#E7298A", "#7570B3", "#1b9e77")

# Common plotting function for violin + beeswarm
plot_violin_beeswarm <- function(data, x, y, y_label, comparisons, pal) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_violin(
      scale = "width", adjust = 0.8,
      draw_quantiles = c(0.25, 0.5, 0.75),
      width = 0.7, alpha = 0.7
    ) +
    geom_beeswarm(cex = 3, alpha = 0.7, size = 1) +
    geom_signif(
      comparisons = comparisons, vjust = 0.2,
      step_increase = 0.1, textsize = 3,
      map_signif_level = FALSE, test = wilcox.test
    ) +
    labs(y = y_label, x = "") +
    scale_fill_manual(values = pal) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(fill = "none")
}

g1 <- plot_violin_beeswarm(data, "type", "ShannonEntropy",
  "Shannon Entropy", comparisons, pal)
g2 <- plot_violin_beeswarm(data, "type", "normalizedShannonEntropy",
  "Normalized Shannon Entropy", comparisons, pal)
g3 <- plot_violin_beeswarm(data, "type", "clonalityIndex",
  "Clonality Index", comparisons, pal)
g4 <- plot_violin_beeswarm(data, "type", "richness",
  "Relative Clonotype Richness", comparisons, pal)

g5 <- plot_violin_beeswarm(data, "type", "reads",
  "Total number of productive TCRs", comparisons, pal) +
  scale_y_log10(limits = c(min(data$reads), NA))

g6 <- plot_violin_beeswarm(data, "type", "diversity",
  "Total number of unique TCRs", comparisons, pal) +
  scale_y_log10(limits = c(min(data$diversity), NA))

# Sample similarity indices
data1 <- read.table("data/TCR/compare_TCR_Statistics.all.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
data1 <- subset(data1, selected == 1)
data1$compare <- factor(data1$compare, levels = c("PvsL", "TvsP", "TvsL"))

comparisons1 <- list(c("PvsL", "TvsP"), c("TvsP", "TvsL"), c("PvsL", "TvsL"))
pal1 <- c("#fdc086", "#d95f02", "#5e797f")

g7 <- plot_violin_beeswarm(data1, "compare", "JaccardIndex",
  "Jaccard Index", comparisons1, pal1)
g8 <- plot_violin_beeswarm(data1, "compare", "HornIndex",
  "Horn Index", comparisons1, pal1)

# Combine all panels
combined <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + plot_layout(nrow = 3)
ggsave("results/TCR_analysis/allsingleT.TCR_BNT.Statistics.pdf",
  combined, width = 7.5, height = 9
)
message("TCR diversity analysis saved.")
