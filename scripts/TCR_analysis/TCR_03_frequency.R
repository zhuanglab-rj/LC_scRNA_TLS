#!/usr/bin/env Rscript
# =============================================================================
# TCR Frequency and Expansion Analysis
#
# Analyzes TCR clone frequency distribution and expansion across sample types:
#   - TCR counts at various frequency thresholds
#   - Expanded TCR fraction comparison
#   - Clonetype proportion analysis
#
# Input:  data/TCR/allsingleT.TCR.percent.matrix.txt
# Output: Multiple intermediate files and results/TCR_analysis/*.pdf
# =============================================================================

source("scripts/00_load_packages.R")

data0 <- read.table("data/TCR/allsingleT.TCR.percent.matrix.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)

# ===========================================================================
# TCR counts at frequency thresholds
# ===========================================================================
sample <- vector()
type <- vector()
percentage <- vector()
number <- vector()

thresholds <- c(0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006)

for (p in thresholds) {
  for (i in 6:ncol(data0)) {
    data <- na.omit(data0[i])
    sample <- c(sample, colnames(data)[1])
    type0 <- strsplit(colnames(data)[1], "_")[[1]][2]
    type0 <- strsplit(type0, "1")[[1]][1]
    type0 <- strsplit(type0, "2")[[1]][1]
    type <- c(type, type0)
    percentage <- c(percentage, p)
    n_above <- sum(data > p)
    number <- c(number, n_above)
  }
}

res1 <- data.frame(sample, type, number, percentage)
write.table(res1, file = "data/TCR/allsingleT.TCR_frequency_threshold.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, dec = "."
)

# Plot frequency threshold analysis
data <- read.table("data/TCR/allsingleT.TCR_frequency_threshold.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
data$percentage <- factor(data$percentage, levels = thresholds)
data$type <- ifelse(data$type == "T", "Tumor",
  ifelse(data$type %in% c("N", "L"), "Lung", "PBMC")
)
data$type <- factor(data$type, levels = c("PBMC", "Lung", "Tumor"))
pal <- c("#E7298A", "#7570B3", "#1b9e77")

p_freq <- ggplot(data, aes(x = percentage, y = (number + 1), fill = type)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0) +
  geom_beeswarm(dodge.width = 0.8, cex = 0.6, size = 1, shape = 21) +
  labs(y = "Number of TCRs above threshold", x = "TCR frequency threshold") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_y_log10()

ggsave("results/TCR_analysis/allsingleT.TCR_frequency_threshold.pdf",
  p_freq, width = 10, height = 5
)

# ===========================================================================
# Expanded TCR fraction analysis
# ===========================================================================
sub0 <- subset(data, percentage == 0)
sub1 <- subset(data, percentage == 0.001)
stopifnot(identical(sub0$sample, sub1$sample))

sub0$expand <- sub1$number
sub0$expand_freq <- sub0$expand / sub0$number

comparisons <- list(c("PBMC", "Lung"), c("PBMC", "Tumor"), c("Lung", "Tumor"))

p_expand <- ggplot(sub0, aes(x = type, y = expand_freq, fill = type)) +
  geom_violin(
    scale = "width", adjust = 0.8,
    draw_quantiles = c(0.25, 0.5, 0.75),
    width = 0.7, alpha = 0.7
  ) +
  geom_beeswarm(cex = 2.5, alpha = 0.7, size = 1) +
  geom_signif(
    comparisons = comparisons, vjust = 0.2,
    step_increase = 0.1, textsize = 3,
    map_signif_level = FALSE, test = wilcox.test
  ) +
  labs(y = "Ratio of expanded/total TCR", x = "") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none")

ggsave("results/TCR_analysis/allsingleT.TCR_expanded_fraction.pdf",
  p_expand, width = 2.5, height = 3
)

# ===========================================================================
# Clonetype proportion analysis
# ===========================================================================
sample <- vector()
patient <- vector()
type <- vector()
N0 <- vector()
N0.001 <- vector()
p0.001 <- vector()
p0.002 <- vector()
p0.0005 <- vector()

for (i in 6:ncol(data0)) {
  data <- na.omit(data0[i])
  sample[i] <- colnames(data)[1]
  patient[i] <- strsplit(sample[i], "_")[[1]][1]
  type[i] <- strsplit(sample[i], "_")[[1]][2]
  p0.001[i] <- sum(data[data > 0.001])
  N0[i] <- sum(data > 0)
  N0.001[i] <- sum(data > 0.001)
  p0.002[i] <- sum(data[data > 0.002])
  p0.0005[i] <- sum(data[data > 0.0005])
}

res <- data.frame(sample, patient, type, N0, N0.001, p0.001, p0.002, p0.0005)[-1:-5, ]
res$type <- ifelse(res$type == "T", "Tumor",
  ifelse(res$type %in% c("N", "L"), "Lung", "PBMC")
)

write.table(res, file = "data/TCR/allsingleT.TCR_clonetype_proportion.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, dec = "."
)

data <- read.table("data/TCR/allsingleT.TCR_clonetype_proportion.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
data$type <- factor(data$type, levels = c("PBMC", "Lung", "Tumor"))

g1 <- ggplot(data, aes(x = type, y = N0.001 / N0, fill = type)) +
  geom_violin(
    scale = "width", adjust = 0.8,
    draw_quantiles = c(0.25, 0.5, 0.75),
    width = 0.7, alpha = 0.7
  ) +
  geom_beeswarm(cex = 2.5, alpha = 0.7, size = 1) +
  geom_signif(
    comparisons = comparisons, vjust = 0.2,
    step_increase = 0.1, textsize = 3,
    map_signif_level = FALSE, test = wilcox.test
  ) +
  labs(y = "ratio of expanded/total TCR", x = "") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none")

g2 <- ggplot(data, aes(x = type, y = p0.001, fill = type)) +
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
  labs(y = "Proportion of TCR repertoire \nrepresented by expanded TCRs", x = "") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  ylim(min(data$p0.001), 1)

data1 <- subset(data, type != "PBMC")
g3 <- ggplot(data1, aes(x = type, y = p0.001, fill = type)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0) +
  geom_beeswarm(cex = 4, size = 1.5, shape = 21) +
  geom_signif(
    comparisons = list(c("Lung", "Tumor")), vjust = 0.2,
    step_increase = 0.1, textsize = 3,
    map_signif_level = FALSE, test = wilcox.test
  ) +
  labs(y = "Proportion of TCR repertoire \nrepresented by expanded TCRs", x = "") +
  scale_fill_manual(values = c("#7570B3", "#1b9e77")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none")

combined <- g1 + g2 + g3 + plot_layout(nrow = 1, width = c(2, 2, 1.5))
ggsave("results/TCR_analysis/allsingleT.TCR_clonetype_ratio&proportion0.001.BNT.pdf",
  combined, width = 7, height = 3
)
message("TCR frequency analysis saved.")
