#!/usr/bin/env Rscript
# =============================================================================
# TCR Clonetype Analysis - Unique vs Shared Clones
#
# Categorizes TCR clones in tumor samples into:
#   - Unique to tumor
#   - Shared with normal lung tissue
#   - Shared with blood (PBMC)
#   - Shared with both normal and blood
#
# Input:  data/TCR/TCR_Statistics.all.txt
#         data/TCR/allsample.TCR.percent.matrix.txt
# Output: data/TCR/allsingleT.TCR.percent.matrix.txt (renamed samples)
#         data/TCR/allsingleT.TCR_tumor_clonetype_unique-share.txt
#         results/TCR_analysis/allsingleT.TCR_tumor_clonetype_unique-share.barplot.pdf
# =============================================================================

source("scripts/00_load_packages.R")

# Load and prepare metadata
meta <- read.table("data/TCR/TCR_Statistics.all.txt",
  header = TRUE, row.names = 1, sep = "\t", as.is = TRUE
)

# Load TCR percent matrix and subset to selected samples
data0 <- read.table("data/TCR/allsample.TCR.percent.matrix.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)
data0 <- data0[, c(colnames(data0)[1:5], rownames(meta))]

# Verify sample order matches
stopifnot(identical(rownames(meta), colnames(data0)[-1:-5]))

# Rename samples using metadata
colnames(data0) <- c(colnames(data0)[1:5], meta$RJLC)
write.table(data0, file = "data/TCR/allsingleT.TCR.percent.matrix.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, dec = "."
)

# Define tumor samples
tumor_samples <- c(
  "RJLC20_T", "RJLC28_T", "RJLC26_T", "RJLC25_T", "RJLC22_T",
  "RJLC21_T", "RJLC24_T", "RJLC267_T", "RJLC268_T", "RJLC269_T",
  "RJLC270_T", "RJLC272_T", "RJLC273_T", "RJLC274_T", "RJLC275_T",
  "RJLC276_T", "RJLC277_T", "RJLC278_T", "RJLC279_T", "RJLC280_T",
  "RJLC281_T", "RJLC283_T"
)

# Classify TCR clones by tissue sharing pattern
sample <- vector()
clonetypeT <- vector()
number <- vector()
Fraction <- vector()
type <- vector()

for (i in seq_along(tumor_samples)) {
  tumor_id <- tumor_samples[i]
  patient <- strsplit(tumor_id, "_")[[1]][1]
  normal_id <- paste0(patient, "_L")
  blood_id <- paste0(patient, "_P")

  data <- data0[, c(tumor_id, normal_id, blood_id)]
  data <- data[!is.na(data[1]) | !is.na(data[2]) | !is.na(data[3]), ]

  # Unique to tumor
  sample <- c(sample, tumor_id)
  clonetypeT <- c(clonetypeT, nrow(data[!is.na(data[1]), ]))
  number <- c(number, nrow(data[!is.na(data[1]) & is.na(data[2]) & is.na(data[3]), ]))
  Fraction <- c(Fraction, sum(data[!is.na(data[1]) & is.na(data[2]) & is.na(data[3]), 1]))
  type <- c(type, "Unique-TCR clone")

  # Shared with normal
  sample <- c(sample, tumor_id)
  clonetypeT <- c(clonetypeT, nrow(data[!is.na(data[1]), ]))
  number <- c(number, nrow(data[!is.na(data[1]) & !is.na(data[2]) & is.na(data[3]), ]))
  Fraction <- c(Fraction, sum(data[!is.na(data[1]) & !is.na(data[2]) & is.na(data[3]), 1]))
  type <- c(type, "Shared-TCR clone with normal")

  # Shared with blood
  sample <- c(sample, tumor_id)
  clonetypeT <- c(clonetypeT, nrow(data[!is.na(data[1]), ]))
  number <- c(number, nrow(data[!is.na(data[1]) & is.na(data[2]) & !is.na(data[3]), ]))
  Fraction <- c(Fraction, sum(data[!is.na(data[1]) & is.na(data[2]) & !is.na(data[3]), 1]))
  type <- c(type, "Shared-TCR clone with blood")

  # Shared with both
  sample <- c(sample, tumor_id)
  clonetypeT <- c(clonetypeT, nrow(data[!is.na(data[1]), ]))
  number <- c(number, nrow(data[!is.na(data[1]) & !is.na(data[2]) & !is.na(data[3]), ]))
  Fraction <- c(Fraction, sum(data[!is.na(data[1]) & !is.na(data[2]) & !is.na(data[3]), 1]))
  type <- c(type, "Shared-TCR clone with normal and blood")
}

res1 <- data.frame(sample, clonetypeT, number, Fraction, type)
write.table(res1, file = "data/TCR/allsingleT.TCR_tumor_clonetype_unique-share.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, dec = "."
)

# Create stacked bar plot
data <- read.table("data/TCR/allsingleT.TCR_tumor_clonetype_unique-share.txt",
  header = TRUE, sep = "\t", as.is = TRUE
)

data$sample <- factor(data$sample,
  levels = c(
    "RJLC269_T", "RJLC277_T", "RJLC281_T", "RJLC272_T", "RJLC268_T",
    "RJLC274_T", "RJLC280_T", "RJLC270_T", "RJLC278_T", "RJLC21_T",
    "RJLC273_T", "RJLC283_T", "RJLC26_T", "RJLC267_T", "RJLC279_T",
    "RJLC28_T", "RJLC20_T", "RJLC275_T", "RJLC276_T", "RJLC24_T",
    "RJLC25_T", "RJLC22_T"
  )
)

data$type <- factor(data$type,
  levels = rev(c(
    "Shared-TCR clone with normal",
    "Shared-TCR clone with normal and blood",
    "Shared-TCR clone with blood",
    "Unique-TCR clone"
  ))
)

myPalette <- rev(c("#5e797f", "#85bb29", "#d95f02", "#1b9e77"))

p <- ggplot(data, aes(y = sample, x = Fraction, fill = type)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_x_continuous(expand = c(0, 0), position = "top") +
  scale_fill_manual(values = myPalette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  ylab("") +
  xlab("Frequence of TCR clones in tumor") +
  guides(fill = guide_legend(reverse = TRUE))

ggsave("results/TCR_analysis/allsingleT.TCR_tumor_clonetype_unique-share.barplot.pdf",
  p, width = 6.5, height = 4
)
message("TCR clonetype analysis saved.")
