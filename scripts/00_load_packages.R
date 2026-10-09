#!/usr/bin/env Rscript
# =============================================================================
# Load all required packages for the LC_scRNA_TLS analysis pipeline
# =============================================================================

# Single-cell RNA-seq analysis
library(Seurat)
library(Matrix)

# Data manipulation
library(dplyr)
library(tidyr)

# Visualization
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggbeeswarm)
library(ggsignif)
library(ggalt)
library(ggalluvial)
library(gridExtra)
library(ggbeeswarm)

# Excel file handling
library(xlsx)

# Heatmap visualization
library(ComplexHeatmap)

# Custom geoms
source("scripts/utils/geom_split_violin.R")

message("All packages loaded successfully.")
