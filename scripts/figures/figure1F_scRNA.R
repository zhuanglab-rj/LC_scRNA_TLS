#!/usr/bin/env Rscript
# =============================================================================
# Figure 1F - scRNA-seq tSNE visualization of lung cancer samples
#
# Panel 1: tSNE plot of cells colored by cluster identity
# Panel 2: tSNE density plots faceted by sample site
#
# Input:  data/res0.4_YKB_dim15.rData (Seurat object)
# Output: results/figures/Figure1F_1.pdf, results/figures/Figure1F_2.pdf
# =============================================================================

source("scripts/00_load_packages.R")

# Load Seurat object
load("data/res0.4_YKB_dim15.rData")

# Panel 1: tSNE plot colored by cluster
p1 <- TSNEPlot(
  object = sc,
  pt.size = 5,
  do.return = TRUE,
  vector.friendly = TRUE,
  plot.title = "RJLC13",
  colors.use = c(
    "#E41A1C", "#FFFF33", "#A65628", "#984EA3",
    "#33A02C", "#386CB0", "#F781BF", "#FF7F00"
  )
)
ggsave("results/figures/Figure1F_1.pdf", p1, width = 5.5, height = 4.3)
message("Figure 1F panel 1 saved.")

# Panel 2: tSNE density plots by sample site
cell_embeddings <- as.data.frame(sc@dr$tsne@cell.embeddings)
cell_embeddings$site <- substr(rownames(cell_embeddings), 1, 4)

p2 <- ggplot(cell_embeddings, aes(x = tSNE_1, y = tSNE_2)) +
  facet_grid(. ~ site) +
  stat_bkde2d(
    bandwidth = c(2, 2),
    aes(fill = ..density..),
    geom = "raster",
    contour = FALSE,
    interpolate = TRUE
  ) +
  geom_point(size = 0.01, color = "white", alpha = 0.1) +
  scale_fill_gradientn(
    colours = colorRampPalette(c(
      "#110E25", "#41266C", "#74287C", "#AE3674",
      "#DF555D", "#ED9469", "#F8D99D"
    ))(14)
  ) +
  theme(panel.background = element_rect(fill = "#110E25"))

ggsave("results/figures/Figure1F_2.pdf", p2, width = 10, height = 5)
message("Figure 1F panel 2 saved to results/figures/Figure1F_2.pdf")
