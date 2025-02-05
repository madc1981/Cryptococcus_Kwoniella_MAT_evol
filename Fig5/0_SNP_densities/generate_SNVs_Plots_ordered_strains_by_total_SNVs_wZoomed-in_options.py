#!/usr/bin/env python3
"""
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2024-12-13
Version: 1.0
Description: This script generates plots for SNV distributions per chromosome and zoomed-in plots for specified regions.
Requirements: pandas, matplotlib, numpy
Usage: python3 Generate_SNVs_Plots_per_chrs.py --merged-snps-file <path> --chr-lengths-file <path> [--regions-file <path>] [--zoom-regions-file <path>] [--region-bar <chromosome:start-end>]
"""

import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib.backends.backend_pdf import PdfPages

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate SNV plots per chromosome and zoomed-in regions.')
parser.add_argument('--merged-snps-file', type=str, required=True, help='Path to the merged SNPs file in BED-like format (Strain, Chromosome, Start, End).')
parser.add_argument('--chr-lengths-file', type=str, required=True, help='Path to the chromosome lengths file (Chromosome, Length).')
parser.add_argument('--regions-file', type=str, default=None, help='Optional path to the regions file (Chromosome, Start, End, Region_Name, Color).')
parser.add_argument('--zoom-regions-file', type=str, default=None, help='Optional file with zoom regions (Chromosome, Start, End).')
parser.add_argument('--region-bar', type=str, default=None, help='Optional region of interest in the format chromosome:start-end.')
args = parser.parse_args()

# File paths
merged_snps_file = args.merged_snps_file
chr_lengths_file = args.chr_lengths_file
regions_file = args.regions_file
zoom_regions_file = args.zoom_regions_file

# Validate optional region input
region_bar = None
if args.region_bar:
    try:
        chrom, positions = args.region_bar.split(':')
        start, end = map(int, positions.split('-'))
        region_bar = {'Chromosome': chrom, 'Start': start, 'End': end}
    except ValueError:
        print("Invalid format for --region-bar. Use 'chromosome:start-end'.")
        exit(1)

# Read data files
merged_df = pd.read_csv(merged_snps_file, sep='\t', header=None, names=['Strain', 'Chromosome', 'Start', 'End'], dtype={'Strain': str}, low_memory=False)
chr_length_df = pd.read_csv(chr_lengths_file, sep='\t', header=None, names=['Chromosome', 'Length'])

# Read regions file if provided
if regions_file:
    regions_df = pd.read_csv(regions_file, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Region_Name', 'Color'])
else:
    regions_df = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Region_Name', 'Color'])

# Read zoom regions file if provided
if zoom_regions_file:
    zoom_regions = pd.read_csv(zoom_regions_file, sep='\t', header=None, names=['Chromosome', 'Start', 'End'])
    if not zoom_regions.empty:
        zoom_regions['Start'] = zoom_regions['Start'].astype(int)
        zoom_regions['End'] = zoom_regions['End'].astype(int)
    else:
        zoom_regions = pd.DataFrame(columns=['Chromosome', 'Start', 'End'])
else:
    zoom_regions = pd.DataFrame(columns=['Chromosome', 'Start', 'End'])

# Count total SNVs for each strain
strain_snv_count = merged_df['Strain'].value_counts().reset_index()
strain_snv_count.columns = ['Strain', 'SNV_Count']
strain_snv_count = strain_snv_count.sort_values(by='SNV_Count', ascending=True)
sorted_strains = strain_snv_count['Strain'].values

# Create colormap for strains
normed_snv_count = (strain_snv_count['SNV_Count'] - strain_snv_count['SNV_Count'].min()) / (strain_snv_count['SNV_Count'].max() - strain_snv_count['SNV_Count'].min())
gradient_colors = plt.cm.BuPu(normed_snv_count)

def plot_snv_distribution(chrom, regions_df, start=None, end=None, suffix=''):
    """Plots SNV distributions for a chromosome or a zoomed-in region."""
    fig, ax = plt.subplots(figsize=(60, 5) if start is None else (30, 5))

    # Filter data for the current chromosome
    chrom_data = merged_df[merged_df['Chromosome'] == chrom]
    chrom_length = chr_length_df[chr_length_df['Chromosome'] == chrom]['Length'].iloc[0]
    regions_for_chrom = regions_df[regions_df['Chromosome'] == chrom] if not regions_df.empty else pd.DataFrame()

    # Apply zoom region if specified
    if start and end:
        chrom_data = chrom_data[(chrom_data['End'] >= start) & (chrom_data['Start'] <= end)]
        chrom_length = end - start  # Update the length for zoomed-in region

    # Scale factors for positions
    scale_factor = 1e6 if chrom_length >= 1e6 else 1e3
    scale_label = 'Mb' if chrom_length >= 1e6 else 'kb'

    # Adjust SNV positions for zoomed-in regions
    scaled_positions = chrom_data['Start'] / scale_factor

    # Group data by strain and plot SNVs
    grouped_data = scaled_positions.groupby(chrom_data['Strain']).apply(list)
    for idx, strain in enumerate(sorted_strains):
        if strain in grouped_data:
            snv_positions = grouped_data[strain]
            ax.vlines(snv_positions, ymin=idx, ymax=idx + 1, color=gradient_colors[np.where(sorted_strains == strain)[0][0]])

    # Add regions of interest as horizontal bars
    for _, region in regions_for_chrom.iterrows():
        region_start = region['Start'] / scale_factor
        region_end = region['End'] / scale_factor
        ax.barh(len(sorted_strains) + 0.3, region_end - region_start, left=region_start, height=0.3,
                color=region['Color'], edgecolor='black', align='center', label=region['Region_Name'])

    # Add optional region bar if specified
    if region_bar and region_bar['Chromosome'] == chrom:
        region_start = region_bar['Start'] / scale_factor
        region_end = region_bar['End'] / scale_factor
        ax.barh(len(sorted_strains) + 0.3, region_end - region_start, left=region_start, height=0.3,
                color='#f1cd1c', edgecolor='black', align='center', label='Region of Interest')

    # Set labels and title
    ax.set_yticks(np.arange(len(sorted_strains)) + 0.5)
    ax.set_yticklabels(sorted_strains, fontsize=12)
    ax.tick_params(axis='x', labelsize=12)

    ax.set_xlabel(f'Length ({scale_label})', fontsize=14)
    ax.set_ylabel('Strain', fontsize=14)
    title = f'{chrom}' if not suffix else f'{chrom} ({start}-{end})'
    ax.set_title(title, fontsize=16)

    # Explicitly set x-axis limits for zoomed-in regions
    if start and end:
        ax.set_xlim(start / scale_factor, end / scale_factor)

    # Add colorbar for SNV counts
    sm = plt.cm.ScalarMappable(cmap='BuPu', norm=plt.Normalize(vmin=strain_snv_count['SNV_Count'].min(), vmax=strain_snv_count['SNV_Count'].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Total SNVs', fontsize=14)

    # Tweak layout and save the figure
    plt.tight_layout()
    plt.savefig(f'SNVs_{chrom}{suffix}.png', bbox_inches='tight')
    plt.savefig(f'SNVs_{chrom}{suffix}.pdf', bbox_inches='tight')
    plt.close(fig)

# Generate plots for each chromosome
for chrom in chr_length_df['Chromosome']:
    plot_snv_distribution(chrom, regions_df)

# Generate zoomed-in plots
for _, row in zoom_regions.iterrows():
    plot_snv_distribution(row['Chromosome'], regions_df, row['Start'], row['End'], suffix=f'_zoom_{row["Start"]}_{row["End"]}')

# Output ordered strains to text file
with open('Ordered_Strains.txt', 'w') as f:
    for strain in sorted_strains:
        f.write(strain + '\n')

