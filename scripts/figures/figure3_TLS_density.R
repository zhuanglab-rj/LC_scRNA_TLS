#!/usr/bin/env Rscript
# =============================================================================
# Figure 3 - TLS density analysis across lesion types
#
# Panel C: TLS density comparison between pGGN and mGGN
# Panel D: TLS density comparison across histological stages (AIS, MIA, IAC)
# Panel E: Paired TLS density comparison (Pre-invasive vs Invasive lesions)
# Panel F: Patient-level TLS density visualization by histology
# Panel G: Alluvial plot showing patient progression across classifications
#
# Input:  data/Fig3C+3D.xlsx, data/Fig3E.txt, data/Fig3F.txt, data/Fig3G.xlsx
# Output: results/figures/Fig3C_3D_plot.pdf, Fig3E_3F_plot.pdf, Fig3G.pdf
# =============================================================================

source("scripts/00_load_packages.R")

# ===========================================================================
# Panels C & D: TLS density by GGO type and histology
# ===========================================================================
data <- read.xlsx("data/Fig3C+3D.xlsx", sheetIndex = 1, row.names = 1)
data$GGO <- factor(data$GGO, levels = c("pGGN", "mGGN"))
data$Histology <- factor(data$Histology, levels = c("AIS", "MIA", "IAC"))

comparisons_ggo <- list(c("pGGN", "mGGN"))
plot_c <- ggplot(data, aes(x = GGO, y = TLS.density, fill = GGO)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.7),
    alpha = 0.4, size = 1.5
  ) +
  labs(y = "TLS density", x = "", tag = "C") +
  scale_fill_manual(values = c("#e6f5c9", "#66c2a5")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.tag = element_text(size = 18, face = "bold")
  ) +
  guides(fill = "none") +
  geom_signif(
    comparisons = comparisons_ggo, vjust = 0.1,
    step_increase = 0.1, textsize = 3,
    map_signif_level = FALSE, test = wilcox.test
  ) +
  ylim(0, 2.35)

comparisons_hist <- list(c("AIS", "MIA"), c("MIA", "IAC"))
plot_d <- ggplot(data, aes(x = Histology, y = TLS.density, fill = Histology)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.9),
    alpha = 0.4, size = 1.5
  ) +
  labs(y = "TLS density", x = "", tag = "D") +
  scale_fill_manual(values = c("#EBCC2A", "#3B9AB2", "#F21A00")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.tag = element_text(size = 18, face = "bold")
  ) +
  guides(fill = "none") +
  geom_signif(
    comparisons = comparisons_hist, vjust = 0.1,
    step_increase = 0.1, textsize = 3,
    map_signif_level = FALSE, test = wilcox.test
  ) +
  ylim(0, 2.35)

plot_cd <- plot_c + plot_d + plot_layout(ncol = 2, width = c(2, 3))
ggsave("results/figures/Fig3C_3D_plot.pdf", plot_cd, width = 5.5, height = 3)
message("Figure 3C-3D saved.")

# ===========================================================================
# Panels E & F: Paired and patient-level TLS density analysis
# ===========================================================================
fig3e <- read.table("data/Fig3E.txt", header = TRUE, sep = "\t", as.is = TRUE)
fig3e_long <- fig3e %>%
  select(Patient, Pre_invasive, Invasive) %>%
  pivot_longer(
    cols = c(Pre_invasive, Invasive),
    names_to = "Lesion", values_to = "TLS_density"
  ) %>%
  mutate(Lesion = factor(Lesion,
    levels = c("Pre_invasive", "Invasive"),
    labels = c("Pre-invasive\nlesions", "Invasive\nlesions")
  ))

# Shapiro-Wilk normality test
normality_pre <- shapiro.test(fig3e$Pre_invasive)$p.value > 0.05
normality_inv <- shapiro.test(fig3e$Invasive)$p.value > 0.05
message("Normality test - Pre-invasive: ", ifelse(normality_pre, "normal", "non-normal"))
message("Normality test - Invasive: ", ifelse(normality_inv, "normal", "non-normal"))

# Wilcoxon paired test
p_val_e <- wilcox.test(fig3e$Pre_invasive, fig3e$Invasive, paired = TRUE)$p.value
p_label_e <- ifelse(p_val_e < 0.001, "P < 0.001", paste0("P = ", signif(p_val_e, 3)))

plot_e <- ggplot(fig3e_long, aes(x = Lesion, y = TLS_density, fill = Lesion)) +
  geom_boxplot(width = 0.35, outlier.shape = NA, color = "black") +
  geom_line(aes(group = Patient), color = "grey60", linewidth = 0.4) +
  geom_point(shape = 21, size = 3.5, color = "black", stroke = 0.6) +
  scale_fill_manual(
    values = c("Pre-invasive\nlesions" = "#CCFFFF", "Invasive\nlesions" = "#1f78b4"),
    guide = "none"
  ) +
  annotate("segment", x = 1, xend = 2, y = 2.15, yend = 2.15, linewidth = 0.5) +
  annotate("text", x = 1.5, y = 2.25, label = p_label_e, size = 4.5) +
  scale_y_continuous(limits = c(-0.05, 2.6), breaks = seq(0, 2.5, 0.5)) +
  labs(x = NULL, y = "TLS density", tag = "E") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.margin = margin(5, 40, 5, 5)
  )

# Panel F: Patient-level dot plot
fig3f <- read.table("data/Fig3F.txt", header = TRUE, sep = "\t", as.is = TRUE)
fig3f <- fig3f %>% mutate(Histology = factor(Histology, levels = c("AIS", "MIA", "IAC")))

patient_order <- fig3f %>%
  filter(Histology == "IAC") %>%
  arrange(TLS_density) %>%
  pull(Patient)
fig3f$Patient <- factor(fig3f$Patient, levels = patient_order)

median_AIS <- median(subset(fig3f, Histology == "AIS")$TLS_density)
median_MIA <- median(subset(fig3f, Histology == "MIA")$TLS_density)
median_IAC <- median(subset(fig3f, Histology == "IAC")$TLS_density)

hist_colors <- c("AIS" = "#EBCC2A", "MIA" = "#3B9AB2", "IAC" = "#F21A00")

plot_f <- ggplot(fig3f, aes(x = TLS_density, y = Patient, group = Patient, fill = Histology)) +
  geom_vline(xintercept = median_AIS, color = "#EBCC2A", linewidth = 0.7) +
  geom_vline(xintercept = median_MIA, color = "#3B9AB2", linewidth = 0.7) +
  geom_vline(xintercept = median_IAC, color = "#F21A00", linewidth = 0.7) +
  geom_line(color = "grey75", linetype = 2) +
  geom_point(aes(fill = Histology), shape = 21, size = 4.5, color = "black", stroke = 0.6) +
  scale_fill_manual(values = hist_colors, name = NULL) +
  labs(x = "TLS density", y = NULL, tag = "F") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black", size = 11),
    axis.line = element_line(color = "black"),
    legend.position = c(0.85, 0.18),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.margin = margin(5, 40, 5, 5)
  )

plot_ef <- plot_e / plot_f + plot_layout(heights = c(1, 1.6))
ggsave("results/figures/Fig3E_3F_plot.pdf", plot_ef, width = 5, height = 8)
message("Figure 3E-3F saved.")

# ===========================================================================
# Panel G: Alluvial plot showing patient classification transitions
# ===========================================================================
df <- read.xlsx("data/Fig3G.xlsx", sheetIndex = 1, row.names = 1)
UCB_lodes <- to_lodes_form(df[, 1:ncol(df)],
  axes = 1:ncol(df),
  id = "Patient"
)

UCB_lodes$stratum <- factor(UCB_lodes$stratum,
  levels = c("pGGN", "mGGN", "AIS", "MIA", "IAC", "0", "0-0.1", "0.1-0.3", ">0.3")
)

mycol <- c(
  "#e6f5c9", "#66c2a5", "#EBCC2A", "#3B9AB2", "#F21A00",
  "#ffffe5", "#fee391", "#fe9929", "#993404"
)

plot_g <- ggplot(UCB_lodes, aes(
  x = x, stratum = stratum, alluvium = Patient,
  fill = stratum, label = stratum
)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 1 / 8) +
  geom_stratum(alpha = 0.9, width = 1 / 6) +
  geom_text(stat = "stratum", size = 2, color = "black") +
  scale_fill_manual(values = mycol) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  ggtitle("")

ggsave("results/figures/Fig3G.pdf", plot_g, width = 4, height = 4)
message("Figure 3G saved.")
