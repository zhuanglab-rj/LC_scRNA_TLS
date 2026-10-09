#!/usr/bin/env Rscript
# =============================================================================
# Figure 4 - TLS subtypes density across histological stages
#
# Compares total TLS, immature TLS (iTLS), and mature TLS (mTLS) density
# across AIS, MIA, and IAC histological stages.
#
# Input:  data/Fig4B+4C+4D.xlsx
# Output: results/figures/Fig4B+4C+4D_plot.pdf
# =============================================================================

source("scripts/00_load_packages.R")

data <- read.xlsx("data/Fig4B+4C+4D.xlsx", sheetIndex = 1, row.names = 1)
data$Histology <- factor(data$Histology, levels = c("AIS", "MIA", "IAC"))

comparisons <- list(c("MIA", "IAC"), c("AIS", "IAC"))
hist_colors <- c("#EBCC2A", "#3B9AB2", "#F21A00")

# Common plotting function for TLS density panels
plot_tls_density <- function(data, y_col, y_label, tag) {
  ggplot(data, aes(x = Histology, y = .data[[y_col]], fill = Histology)) +
    geom_boxplot(width = 0.7, outlier.shape = NA) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.9),
      alpha = 0.4, size = 1.5
    ) +
    labs(y = y_label, x = "", tag = tag) +
    scale_fill_manual(values = hist_colors) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.tag = element_text(size = 18, face = "bold")
    ) +
    guides(fill = "none") +
    geom_signif(
      comparisons = comparisons, vjust = 0.1,
      step_increase = 0.1, textsize = 3,
      map_signif_level = FALSE, test = wilcox.test
    ) +
    ylim(0, 2.3)
}

plot_b <- plot_tls_density(data, "TLS.density", "TLS density", "B")
plot_c <- plot_tls_density(data, "iTLS.density", "iTLS density", "C")
plot_d <- plot_tls_density(data, "mTLS.density", "mTLS density", "D")

combined <- plot_b + plot_c + plot_d + plot_layout(ncol = 3)
ggsave("results/figures/Fig4B+4C+4D_plot.pdf", combined, width = 9, height = 3)
message("Figure 4B-4D saved to results/figures/Fig4B+4C+4D_plot.pdf")
