#!/usr/bin/env Rscript
# =============================================================================
# TCR Clone Heatmap Analysis
#
# Generates heatmaps showing TCR clone frequency distribution across samples,
# grouped by patient and tissue type.
#
# Input:  data/TCR/TCR_Statistics.all.txt
#         data/TCR/allsingleT.TCR.percent.matrix.txt
# Output: results/TCR_analysis/allsingleT.TCR.percent0.001.heatmap.pdf
# =============================================================================

source("scripts/00_load_packages.R")

# Load metadata and TCR frequency matrix
pdata <- read.table("data/TCR/TCR_Statistics.all.txt",
  header = TRUE, row.names = 1, sep = "\t", as.is = TRUE
)
rownames(pdata) <- pdata$RJLC

data0 <- read.table("data/TCR/allsingleT.TCR.percent.matrix.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)

stopifnot(identical(rownames(pdata), colnames(data0)[-1:-5]))

# Filter clones with frequency >= 0.001
data0[is.na(data0)] <- 0
data0$max <- apply(data0[, 6:ncol(data0)], 1, max)
data <- subset(data0, max >= 0.001)
score <- data[, c(-1:-5, -ncol(data))]
score <- as.matrix(score)

# Annotation for sample types
annotation_col <- data.frame(type = pdata$type)
rownames(annotation_col) <- colnames(score)
ann_colors <- list(
  type = c("PBMC" = "#E7298A", "Lung" = "#7570B3", "Tumor" = "#1b9e77")
)

# Heatmap 1: ordered by patient
patient_order <- c(
  "RJLC280", "RJLC26", "RJLC276", "RJLC278", "RJLC268", "RJLC267",
  "RJLC25", "RJLC270", "RJLC281", "RJLC20", "RJLC283", "RJLC274",
  "RJLC22", "RJLC277", "RJLC275", "RJLC28", "RJLC272", "RJLC273",
  "RJLC24", "RJLC269", "RJLC21", "RJLC279"
)
pdata$patient <- factor(pdata$patient, levels = patient_order)
pdata <- pdata[order(pdata$patient), ]
score <- score[, match(rownames(pdata), colnames(score))]

annotation_col <- data.frame(type = pdata$type, patient = pdata$patient)
rownames(annotation_col) <- colnames(score)

ggsave("results/TCR_analysis/allsingleT.TCR.percent0.001.heatmap.pdf",
  plot = {
    pdf("results/TCR_analysis/allsingleT.TCR.percent0.001.heatmap.pdf",
      width = 12, height = 10
    )
    ComplexHeatmap::pheatmap(score,
      color = colorRampPalette(c(
        "#08519c", "#08519c", "#2171b5",
        "white", "#ef3b2c", "#cb181d", "#cb181d"
      ))(75),
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      cluster_cols = FALSE,
      cellwidth = 10,
      cellheight = 0.1,
      fontsize_col = 10,
      fontsize_row = 0.1,
      scale = "row",
      column_split = annotation_col$patient
    )
    dev.off()
  }
)
message("TCR heatmap saved.")
