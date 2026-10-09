#!/usr/bin/env Rscript
# =============================================================================
# Figure 2B - Resident vs Circulating T cell characterization
#
# Eight-panel figure showing tSNE plots colored by:
#   - Resident / Circulating T cell labels
#   - Resident markers: CCR5, CXCR6, ITGA1
#   - Circulating markers: CCR7, SELL, S1PR1
#
# Input:  data/dim35_Resident_Circulating.rData (Seurat object)
# Output: results/figures/Figure2B.pdf
# =============================================================================

source("scripts/00_load_packages.R")

# Color palettes
palette <- colorRampPalette(c(
  "#00004B", "#004B96", "#32C8FF", "#69C864",
  "#FFFF00", "orange", "#BF0000"
))(100)
greens <- colorRampPalette(c(
  "white", "#bae4b3", "#74cd76",
  "#31a354", "#006d2c", "#006d2c"
))(100)
reds <- colorRampPalette(c("#ffff99", "red"))(100)

# Load Seurat object with meta data
load("data/dim35_Resident_Circulating.rData")
meta <- sc@meta.data

# Common theme for all panels
base_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_fixed(ratio = 1)

# Panel 1: Resident score
g1 <- ggplot(meta[order(meta$res), ], aes(x = tSNE_1, y = tSNE_2, color = res)) +
  labs(x = "", y = "", title = "Resident") +
  geom_point(size = 1) +
  scale_color_gradientn(colours = palette) + base_theme

# Panel 2: CCR5 (resident marker)
g2 <- ggplot(meta[order(meta$CCR5), ], aes(x = tSNE_1, y = tSNE_2, fill = CCR5)) +
  labs(x = "", y = "", title = "CCR5") +
  geom_point(size = 1, shape = 21, color = "lightgrey", stroke = 0.2) +
  scale_fill_gradientn(colours = greens) + base_theme

# Panel 3: CXCR6 (resident marker)
g3 <- ggplot(meta[order(meta$CXCR6), ], aes(x = tSNE_1, y = tSNE_2, fill = CXCR6)) +
  labs(x = "", y = "", title = "CXCR6") +
  geom_point(size = 1, shape = 21, color = "lightgrey", stroke = 0.2) +
  scale_fill_gradientn(colours = greens) + base_theme

# Panel 4: ITGA1 (resident marker)
g4 <- ggplot(meta[order(meta$ITGA1), ], aes(x = tSNE_1, y = tSNE_2, fill = ITGA1)) +
  labs(x = "", y = "", title = "ITGA1") +
  geom_point(size = 1, shape = 21, color = "lightgrey", stroke = 0.2) +
  scale_fill_gradientn(colours = greens) + base_theme

# Panel 5: Circulating score
g5 <- ggplot(meta[order(meta$cir), ], aes(x = tSNE_1, y = tSNE_2, color = cir)) +
  labs(x = "", y = "", title = "Circulating") +
  geom_point(size = 1) +
  scale_color_gradientn(colours = palette) + base_theme

# Panel 6: CCR7 (circulating marker)
g6 <- ggplot(meta[order(meta$CCR7), ], aes(x = tSNE_1, y = tSNE_2, fill = CCR7)) +
  labs(x = "", y = "", title = "CCR7") +
  geom_point(size = 1, shape = 21, color = "lightgrey", stroke = 0.2) +
  scale_fill_gradientn(colours = reds) + base_theme

# Panel 7: SELL (circulating marker)
g7 <- ggplot(meta[order(meta$SELL), ], aes(x = tSNE_1, y = tSNE_2, fill = SELL)) +
  labs(x = "", y = "", title = "SELL") +
  geom_point(size = 1, shape = 21, color = "lightgrey", stroke = 0.2) +
  scale_fill_gradientn(colours = reds) + base_theme

# Panel 8: S1PR1 (circulating marker)
g8 <- ggplot(meta[order(meta$S1PR1), ], aes(x = tSNE_1, y = tSNE_2, fill = S1PR1)) +
  labs(x = "", y = "", title = "S1PR1") +
  geom_point(size = 1, shape = 21, color = "lightgrey", stroke = 0.2) +
  scale_fill_gradientn(colours = reds) + base_theme

# Combine panels
combined <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + plot_layout(ncol = 4)
ggsave("results/figures/Figure2B.pdf", combined, width = 20, height = 12)
message("Figure 2B saved to results/figures/Figure2B.pdf")
