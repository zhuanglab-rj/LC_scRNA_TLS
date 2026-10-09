#!/usr/bin/env Rscript
# =============================================================================
# TCR Viral Specificity Analysis
#
# Analyzes viral vs non-viral TCR proportions across sample types:
#   - Viral/Non-viral TCR ratio by tissue type
#   - Comparison between expanded and non-expanded TCRs
#   - Analysis by tissue sharing pattern
#
# Prerequisite: perl TCR.clonetype.percent0.001.pl
# Input:  Multiple TCR clonetype statistics files
# Output: results/TCR_analysis/allsingleT.TCR_clonetype_TLPviral_statistics.pdf
# =============================================================================

source("scripts/00_load_packages.R")

# Load metadata
meta <- read.table("data/TCR/TCR_Statistics.all.txt",
  header = TRUE, row.names = 1, sep = "\t", as.is = TRUE
)

# ===========================================================================
# Panel 1: Viral/Non-viral ratio by tissue type
# ===========================================================================
df1 <- read.table("data/TCR/TCR_clonetype_TNBviral_statistics2.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
df2 <- read.table("data/TCR/TCR_clonetype_TLPviral_statistics2.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
data0 <- rbind(df1, df2)
data0 <- data0[data0$sample %in% rownames(meta), ]

write.table(data0, file = "data/TCR/allsingleT.TCR_clonetype_TLPviral_statistics.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, dec = "."
)

data0$type <- factor(data0$type, levels = c("PBMC", "Lung", "Tumor"))
pal <- c("#E7298A", "#7570B3", "#1b9e77")
comparisons <- list(c("PBMC", "Lung"), c("PBMC", "Tumor"), c("Lung", "Tumor"))

g1 <- ggplot(data0, aes(x = type, y = viralNoRatio, fill = type)) +
  geom_violin(
    scale = "width", adjust = 1,
    draw_quantiles = c(0.25, 0.5, 0.75),
    width = 0.7, alpha = 0.8
  ) +
  geom_beeswarm(cex = 4, alpha = 0.7, size = 1) +
  geom_signif(
    comparisons = comparisons, vjust = 0.2,
    step_increase = 0.1, textsize = 3,
    map_signif_level = FALSE, test = wilcox.test
  ) +
  labs(y = "Viral/Non-viral TCRs", x = "") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  scale_y_continuous(trans = "log2")

# ===========================================================================
# Panel 2: Viral ratio in expanded vs non-expanded TCRs
# ===========================================================================
df1 <- read.table("data/TCR/TCR_clonetype_TNBviral_ExpandedTCR2.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
df2 <- read.table("data/TCR/TCR_clonetype_TLPviral_ExpandedTCR2.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
data <- rbind(df1, df2)
data <- data[data$sample %in% rownames(meta), ]

write.table(data, file = "data/TCR/allsingleT.TCR_clonetype_TLPviral_ExpandedTCR.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, dec = "."
)

data$expanded <- factor(data$expanded, levels = c("Non-expanded", "Expanded"))
data$type <- factor(data$type, levels = c("PBMC", "Lung", "Tumor"))

g2 <- ggplot(data, aes(x = expanded, y = viralNoRatio, fill = expanded)) +
  geom_violin(
    scale = "width", adjust = 1,
    draw_quantiles = c(0.25, 0.5, 0.75),
    width = 0.7, alpha = 0.8
  ) +
  geom_beeswarm(cex = 4, alpha = 0.7, size = 1) +
  geom_signif(
    comparisons = list(c("Non-expanded", "Expanded")),
    vjust = 0.2, step_increase = 0.1, textsize = 3,
    map_signif_level = FALSE, test = wilcox.test
  ) +
  labs(y = "Viral/Non-viral TCRs", x = "") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = "none") +
  scale_y_continuous(trans = "log2") +
  facet_grid(. ~ type)

# ===========================================================================
# Panel 3: Viral ratio by tissue sharing pattern
# ===========================================================================
df1 <- read.table("data/TCR/TCR_clonetype_TNBviral_statistics3.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
df2 <- read.table("data/TCR/TCR_clonetype_TLPviral_statistics3.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
data1 <- rbind(df1, df2)
data1 <- data1[data1$patient %in% rownames(meta), ]

write.table(data1, file = "data/TCR/allsingleT.TCR_clonetype_TLPviral_shareTCR.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, dec = "."
)

data1$type <- factor(data1$type,
  levels = c(
    "PBMC unique", "Lung unique", "Tumor unique",
    "PBMC/Lung", "PBMC/Tumor", "Lung/Tumor", "PBMC/Lung/Tumor"
  )
)
pal1 <- c(
  "#E7298A", "#7570B3", "#1b9e77",
  "#fdc086", "#d95f02", "#5e797f", "#85bb29"
)

g3 <- ggplot(data1, aes(x = type, y = viralNoRatio, fill = type)) +
  geom_violin(
    scale = "width", adjust = 1,
    draw_quantiles = c(0.25, 0.5, 0.75),
    width = 0.7, alpha = 0.8
  ) +
  geom_beeswarm(cex = 2, alpha = 0.7, size = 1) +
  labs(y = "Viral/Non-viral TCRs", x = "") +
  scale_fill_manual(values = pal1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0, size = 8)
  ) +
  guides(fill = "none") +
  scale_y_continuous(trans = "log2")

# Combine panels
combined <- g1 + g2 + g3 + plot_layout(nrow = 1, width = c(1, 2.25, 2))
ggsave("results/TCR_analysis/allsingleT.TCR_clonetype_TLPviral_statistics.pdf",
  combined, width = 10.5, height = 4
)
message("TCR viral specificity analysis saved.")
