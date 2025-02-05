#!/usr/bin/env Rscript

# Metadata
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-12-03
# Version: 1.1
# Description: 
#   Analyzes gene presence frequency in MAT loci with the following assumptions:
#     - Pseudogenes are treated as present.
#     - Unclear cases ("?") are treated as absent.
#     - SXI1 and SXI2 are combined into a single column, SXI_combined.
#   Generates two bar plots:
#     1. Gene presence frequency in the preserved input order (SXI_combined included).
#     2. Gene presence frequency sorted by frequency (SXI_combined included).
# Requirements: R packages - dplyr, ggplot2, reshape2
# Usage:
#   1. Place the input file "MAT_loci_genes_presence_absence_matrix_to_plot.tsv" in the working directory.
#   2. Set the working directory using `setwd()` if needed.
#   3. Run this script in an R environment.
# Input:
#   - A tab-delimited file "MAT_loci_genes_presence_absence_matrix_to_plot.tsv" containing a gene presence/absence matrix.
#     The file should include a final row for "Essentiality" classification.
# Output:
#   - Two PDF files:
#       - "{date}_gene_presence_frequency_preserved_order_SXI_combined.pdf": Frequency plot in the original order.
#       - "{date}_gene_presence_frequency_sorted_SXI_combined.pdf": Frequency plot sorted by frequency.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load required libraries
library(dplyr)
library(reshape2)
library(ggplot2)

# Set working directory (if applicable)
#setwd(" ")

# Define today's date in YYYY-MM-DD format
today_date <- format(Sys.Date(), "%Y-%m-%d")

# Load the dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read.delim("MAT_loci_genes_presence_absence_matrix_to_plot.tsv", header = TRUE, check.names = FALSE)

# Clean data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Remove the "X" column (if exists) and trailing unnamed columns
data <- data[, !grepl("^X|Unnamed", colnames(data))]

# Subset matrix data (exclude essentiality row)
matrix_data <- data[-nrow(data), ]  # Exclude the last row (essentiality)

# Fix empty column names
colnames(matrix_data)[colnames(matrix_data) == ""] <- "Unnamed_Column"  # Rename empty column

# Subset presence/absence data
presence_data <- matrix_data[, setdiff(colnames(matrix_data), c("MAT", "Strains", "Unnamed_Column"))]

# Add SXI_combined column: Presence if either SXI1 or SXI2 is present
presence_data$SXI_combined <- apply(matrix_data[, c("SXI1", "SXI2")], 1, function(x) {
  if ("1" %in% x || "pseudo" %in% x) "1" else "0"
})

# Remove SXI1 and SXI2 columns, keeping only SXI_combined
presence_data <- presence_data %>% select(-SXI1, -SXI2)

# Recode data: Treat 'pseudo' as '1' (present) and '?' as '0' (absent)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
recode_data <- presence_data %>% 
  mutate(across(everything(), ~ case_when(
    . == "1" ~ "1",        # Gene present
    . == "pseudo" ~ "1",   # Pseudogene treated as present
    . == "?" ~ "0",        # Unclear treated as absent
    TRUE ~ .               # Keep other values as they are
  )))

# Calculate gene presence frequency
# Summarize presence frequency for each gene (as a proportion of strains)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
presence_summary <- recode_data %>% 
  summarise(across(everything(), ~ sum(. == "1") / nrow(recode_data))) %>% 
  melt(variable.name = "Gene", value.name = "Frequency") # Reshape to long format

# Adjust gene order to match input dataset
# Maintain the original gene order and ensure SXI_combined is included
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
original_gene_order <- colnames(presence_data)  # Original gene order (without SXI1 and SXI2)
adjusted_gene_order <- unique(c("SXI_combined", original_gene_order[!original_gene_order %in% c("SXI1", "SXI2")]))
presence_summary$Gene <- factor(presence_summary$Gene, levels = adjusted_gene_order)

# Define custom colors for HD- and P/R-associated genes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hd_genes <- c("SXI_combined", "RPL22", "CAP1", "SPO14")  # Include SXI_combined as an HD gene
pr_color <- "#58b1a9"  # Color for P/R-associated genes
hd_color <- "#d7b366"  # Color for HD-associated genes

# Assign colors based on gene type
presence_summary <- presence_summary %>%
  mutate(Color = ifelse(Gene %in% hd_genes, hd_color, pr_color))

# Shared theme helper for common formatting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shared_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),  # Remove vertical grid lines
      panel.grid.major.x = element_blank(),  # Remove horizontal grid lines
      panel.grid.minor = element_blank(),    # Remove minor grid lines
      axis.text.y = element_text(size = 9, margin = margin(r = 5)),  # Tick labels for frequency
      axis.title.x = element_text(size = 12),  # Title for x-axis
      axis.title.y = element_text(size = 12),  # Title for y-axis
      plot.title = element_text(size = 14, hjust = 0.5),  # Centered plot title
      axis.ticks.y = element_line(size = 0.5), # Add ticks for y-axis
      axis.ticks.length.y = unit(0.2, "cm"),   # Shorten y-axis tick marks
      axis.line.y = element_line(size = 0.5)   # Add a y-axis bar
    )
}

# Plot 1: Preserved Gene Order
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_preserved_order <- ggplot(presence_summary, aes(x = Gene, y = Frequency, fill = Color)) +
  geom_col() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray75", size = 0.5) +
  labs(
    title = "Gene presence frequency in MAT loci (preserved input order, SXI1/2 combined)",
    x = "",
    y = "Presence frequency"
  ) +
  shared_theme() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, margin = margin(t = -5))) +  # Bottom alignment
  scale_fill_identity()

# Save Plot 1
ggsave(paste0(today_date, "_gene_presence_frequency_preserved_order_SXI_combined.pdf"), plot = plot_preserved_order, width = 10, height = 8, dpi = 300)

# Plot 2: Sorted by Frequency
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_sorted_frequency <- ggplot(presence_summary, aes(x = reorder(Gene, -Frequency), y = Frequency, fill = Color)) + 
  geom_col() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray75", size = 0.5) +
  labs(
    title = "Gene presence frequency in MAT loci (sorted by frequency, SXI1/2 combined)",
    x = "",
    y = "Presence frequency"
  ) +
  shared_theme() +
  theme(axis.text.x.top = element_text(size = 10, angle = 90, hjust = 0, margin = margin(t = -5))) +  # Top alignment
  scale_fill_identity() +
  scale_y_reverse(breaks = seq(0, 1, by = 0.1), limits = c(1, 0)) +
  scale_x_discrete(position = "top")

# Save Plot 2
ggsave(paste0(today_date, "_gene_presence_frequency_sorted_SXI_combined.pdf"), plot = plot_sorted_frequency, width = 10, height = 8, dpi = 300)
