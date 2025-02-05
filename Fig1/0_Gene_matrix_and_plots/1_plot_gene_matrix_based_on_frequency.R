#!/usr/bin/env Rscript

# Metadata
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2024-12-05
# Version: 1.1
# Description: 
#   Generates two heatmaps:
#   1. A heatmap with genes ordered by frequency and SXI_combined (a combination of SXI1 and SXI2).
#   2. A heatmap with genes ordered by frequency, showing SXI1 and SXI2 as separate columns 
#      in the position of SXI_combined.
# Requirements: R packages - dplyr, ComplexHeatmap, circlize
# Usage:
#   1. Place the input file "MAT_loci_genes_presence_absence_matrix_to_plot.tsv" in the working directory.
#   2. Set the working directory using `setwd()`.
#   3. Run this script in an R environment.
# Input:
#   - A tab-delimited file "MAT_loci_genes_presence_absence_matrix_to_plot.tsv" containing the gene presence/absence matrix.
# Output:
#   - Two PDF files:
#       - "{date}_gene_presence_with_SXI_combined.pdf"
#       - "{date}_gene_presence_with_SXI1_SXI2.pdf"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set working directory
# setwd(" ")

# Load required libraries
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Load the dataset
data <- read.delim("MAT_loci_genes_presence_absence_matrix_to_plot.tsv", header = TRUE, check.names = FALSE)

# Data cleaning and preparation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Remove metadata and fix column names
data <- data[, !grepl("^X|Unnamed", colnames(data))]  # Remove unwanted columns
essentiality <- tail(data, 1)  # Extract essentiality row (last row)
matrix_data <- data[-nrow(data), ]  # Exclude essentiality row from the matrix
colnames(matrix_data)[colnames(matrix_data) == ""] <- "Unnamed_Column"  # Fix empty column names

# Subset presence/absence data and add SXI_combined
presence_data <- matrix_data[, setdiff(colnames(matrix_data), c("MAT", "Strains", "Unnamed_Column"))]
presence_data$SXI_combined <- apply(matrix_data[, c("SXI1", "SXI2")], 1, function(x) {
  if ("1" %in% x || "pseudo" %in% x) "1" else "0"
})

# Recode data for numeric analysis
presence_data_colored <- presence_data  # Retain original for coloring
recode_data <- presence_data %>% mutate(across(everything(), ~ case_when(
  . == "1" ~ 1,        # Gene present
  . == "pseudo" ~ 1,   # Pseudogene treated as present
  . == "?" ~ 0,        # Unclear treated as absent
  TRUE ~ 0             # All others treated as absent
)))

# Calculate gene frequency and order
gene_presence <- colSums(recode_data)
ordered_genes <- names(sort(gene_presence, decreasing = TRUE))  # Order genes by presence frequency

# Define shared visual parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
heatmap_colors <- c("1" = "#9a8479", "pseudo" = "#cbb39a", "0" = "#e6e7e8", "?" = "#a6a8ab")
essentiality_colors <- c("ESS" = "#ed1c24", "NESS" = "#1744ff", "UNK" = "#a6a8ab")
genus_colors <- c("Cryptococcus" = "#f8beb1", "Kwoniella" = "#bae3f5")

# Row Annotation: Genus-based coloring
row_annotation <- rowAnnotation(
  Genus = ifelse(grepl("^Cryptococcus_", matrix_data$Strains), "Cryptococcus", "Kwoniella"),
  col = list(Genus = genus_colors),
  annotation_legend_param = list(
    Genus = list(labels_gp = gpar(fontface = "italic"), title = NULL)
  ),
  show_annotation_name = FALSE
)

# Heatmap 1: SXI_combined
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ordered_genes_combined <- ordered_genes[!ordered_genes %in% c("SXI1", "SXI2")]  # Exclude SXI1 and SXI2
essentiality$SXI_combined <- ifelse(essentiality$SXI1 == "ESS" | essentiality$SXI2 == "ESS", "ESS", "NESS")
essentiality_combined <- essentiality %>% select(-SXI1, -SXI2)  # Remove SXI1 and SXI2 from essentiality
essentiality_combined <- unlist(essentiality_combined[, ordered_genes_combined, drop = TRUE])  # Reorder essentiality

recode_matrix_combined <- as.matrix(recode_data[, ordered_genes_combined])
rownames(recode_matrix_combined) <- matrix_data$Strains
colored_matrix_combined <- as.matrix(presence_data_colored[, ordered_genes_combined])
rownames(colored_matrix_combined) <- matrix_data$Strains

# Column Annotation for combined plot
column_annotation_combined <- HeatmapAnnotation(
  Essentiality = essentiality_combined,
  col = list(Essentiality = essentiality_colors),
  annotation_legend_param = list(
    labels = c("Essential", "Non-essential", "Unknown"),
    at = c("ESS", "NESS", "UNK"),
    title = NULL
  ),
  which = "column",
  show_annotation_name = FALSE
)

# Produce combined heatmap
ht_no_clustering_combined <- Heatmap(
  colored_matrix_combined,
  name = " ",
  col = heatmap_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_labels = rownames(colored_matrix_combined),
  show_row_names = TRUE,
  show_column_names = TRUE,
  bottom_annotation = column_annotation_combined,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8, fontface = "italic"),
  height = unit(15, "cm"),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "white", lwd = 1))
  }
) + row_annotation

# Save combined heatmap
today_date <- format(Sys.Date(), "%Y-%m-%d")
pdf(paste0(today_date, "_gene_presence_with_SXI1-2_combined.pdf"), width = 12, height = 15)
draw(ht_no_clustering_combined, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# Heatmap 2: SXI1 and SXI2 separated
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sxi_combined_position <- which(ordered_genes == "SXI_combined")
if (length(sxi_combined_position) > 0) {
  ordered_genes_separated <- c(
    ordered_genes[1:(sxi_combined_position - 1)],  
    "SXI1", "SXI2",                               
    ordered_genes[(sxi_combined_position + 1):length(ordered_genes)]
  )
  ordered_genes_separated <- unique(ordered_genes_separated)  # Ensure no duplicates
}

recode_matrix_separated <- as.matrix(recode_data[, ordered_genes_separated])
rownames(recode_matrix_separated) <- matrix_data$Strains
colored_matrix_separated <- as.matrix(presence_data_colored[, ordered_genes_separated])
rownames(colored_matrix_separated) <- matrix_data$Strains
essentiality_separated <- unlist(essentiality[, ordered_genes_separated, drop = TRUE])

column_annotation_separated <- HeatmapAnnotation(
  Essentiality = essentiality_separated,
  col = list(Essentiality = essentiality_colors),
  annotation_legend_param = list(
    labels = c("Essential", "Non-essential", "Unknown"),
    at = c("ESS", "NESS", "UNK"),
    title = NULL
  ),
  which = "column",
  show_annotation_name = FALSE
)

# Produce separated heatmap
ht_no_clustering_separated <- Heatmap(
  colored_matrix_separated,
  name = " ",
  col = heatmap_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_labels = rownames(colored_matrix_separated),
  show_row_names = TRUE,
  show_column_names = TRUE,
  bottom_annotation = column_annotation_separated,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8, fontface = "italic"),
  height = unit(15, "cm"),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "white", lwd = 1))
  }
) + row_annotation

# Save separated heatmap
pdf(paste0(today_date, "_gene_presence_with_separate_SXI1_SXI2.pdf"), width = 12, height = 15)
draw(ht_no_clustering_separated, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()