#!/usr/bin/env Rscript
# =============================================================================
# Figure 1E - IHC staining quantification
#
# Compares immunohistochemistry (IHC) marker density between tumor and normal
# tissue using split violin plots. Markers include CD3, CD4, CD8, CD20, CD68,
# CD56, MCT, PD1, and PD-L1.
#
# Input:  data/1E.txt (tab-delimited with columns: gene, value2, type)
# Output: results/figures/Figure1E.pdf
# =============================================================================

source("scripts/00_load_packages.R")

# Read IHC quantification data
data <- read.table("data/1E.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)

# Set factor levels for consistent ordering
data$gene <- factor(data$gene,
  levels = c("CD3", "CD4", "CD8", "CD20", "CD68", "CD56", "MCT", "PD1", "PD-L1")
)

# Create split violin plot comparing tumor vs normal tissue
p <- ggplot(data, aes(x = gene, y = log2(value2), fill = type)) +
  geom_split_violin(
    draw_quantiles = c(0.25, 0.5, 0.75),
    trim = FALSE,
    linetype = "solid",
    color = "black",
    alpha = 0.7,
    size = 0.5,
    na.rm = TRUE,
    position = "identity"
  ) +
  ylab("IHC [log(1+density)]") +
  xlab("") +
  scale_fill_manual(values = c("#7570B3", "#1b9e77")) +
  theme_set(theme_set(theme_classic(base_size = 20))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  guides(fill = guide_legend(title = NULL)) +
  theme(
    legend.background = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0)
  )

# Save output
ggsave("results/figures/Figure1E.pdf", p, width = 6, height = 4)
message("Figure 1E saved to results/figures/Figure1E.pdf")
