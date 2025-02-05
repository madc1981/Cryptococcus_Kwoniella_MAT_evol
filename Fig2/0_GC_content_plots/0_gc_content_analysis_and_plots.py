#!/usr/bin/env python3
"""
Author:         Marco A. Coelho @ Heitman lab, Duke University
Date:           2024-07-03
Version:        1.3
Description:    This script calculates the GC content and GC content deviation from the genome mean GC% using a sliding window approach.
                It generates consolidated plots for each genome, a single raw data file per genome, and optionally zoomed-in views
                with customizable sliding window parameters. It can also highlight specific regions defined in a BED file.

Requirements:   biopython, matplotlib, pandas, argparse, glob
Usage:          python gc_content_analysis.py --window_size 1000 --step_size 1000
                python gc_content_analysis.py --coordinates_file coordinates.txt --zoom_window_size 100 --zoom_step_size 10 --bed_file regions.bed --num_cpus 4 --output_dir'my_output'

Inputs:
    1. Genome FASTA files (*.fa):
        - The script processes all FASTA files with the .fa extension in the current directory.
        - Each FASTA file should contain the genome sequence(s) in standard FASTA format.

    2. Coordinates file (optional, specified with `--coordinates_file`):
        - A tab-delimited file with the following columns:
          [Strain_Species_ID, Contig_ID, Start, End]
        - Used for generating zoomed-in plots of specific genome regions.

    3. BED file (optional, specified with `--bed_file`):
        - A tab-delimited file containing regions to highlight in the plots.
        - Expected columns:
          [Strain/Species_ID, Contig_ID, Start, End, Region_ID, Color_Code]

Outputs:
    1. Consolidated plots of GC content and GC content deviation for each genome.
    2. Individual plots for zoomed-in regions (if coordinates file is provided).
    3. CSV files containing GC content and deviation data for each genome.
    4. All output files are saved in the specified output directory (default: 'output').
"""


import os
import glob
import argparse
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

# Set PDF font type
plt.rcParams['pdf.fonttype'] = 42

def calculate_gc_content(sequence):
    """Calculate GC content of a given sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence) * 100

def sliding_window_gc(seq, window_size, step_size, start=0):
    """Calculate GC content in sliding windows across the sequence."""
    gc_content = []
    positions = []
    for i in range(0, len(seq) - window_size + 1, step_size):
        window_seq = seq[i:i + window_size]
        gc_content.append(calculate_gc_content(window_seq))
        positions.append((start + i + window_size // 2) / 1e6)  # Convert positions to Mb
    return gc_content, positions

def plot_consolidated_gc(genome_data, output_file, y_label, bed_data=None, strain_species_id=None, y_limits=None, figsize=(15, 8)):
    """Plot consolidated GC content and save the plot."""
    fig, axes = plt.subplots(len(genome_data), 1, figsize=figsize, sharex=True)  # Use figsize parameter
    if len(genome_data) == 1:
        axes = [axes]

    # Determine global Y-axis limits if not provided
    if y_limits is None:
        max_gc = max(max(values) for positions, values in genome_data.values())
        min_gc = min(min(values) for positions, values in genome_data.values())
    else:
        min_gc, max_gc = y_limits

    for ax, (contig_id, (positions, values)) in zip(axes, genome_data.items()):
        # Plot the data
        ax.plot(positions, values, label=contig_id, linewidth=0.5)

        # Set contig labels on the left side
        ax.set_ylabel(contig_id, rotation=0, labelpad=20, verticalalignment='center')

        # Set consistent Y-axis limits
        ax.set_ylim(min_gc, max_gc)

        # Move only the ticks to the right side
        ax.yaxis.tick_right()

        # Highlight regions from BED file if provided
        if bed_data is not None and strain_species_id in bed_data and contig_id in bed_data[strain_species_id]:
            for region in bed_data[strain_species_id][contig_id]:
                ax.axvspan(region[0] / 1e6, region[1] / 1e6, color=region[3], alpha=0.5)  # Use specific color code

    # Set X-axis label and figure title
    plt.xlabel("Genomic position (Mb)")
    fig.suptitle(f"{strain_species_id} - GC Content")

    # Add Y-axis label on the right side of the figure closer to the plot area
    fig.text(0.95, 0.5, y_label, va='center', ha='center', rotation='vertical')

    # Save plot as PDF
    plt.savefig(output_file)
    plt.close()

def plot_consolidated_gc_deviation(genome_data, output_file, step_size, bed_data=None, strain_species_id=None, y_limits=None, figsize=(15, 8)):
    """Plot consolidated GC content deviation with colored positive and negative values and save the plot."""
    fig, axes = plt.subplots(len(genome_data), 1, figsize=figsize, sharex=True)  # Use figsize parameter
    if len(genome_data) == 1:
        axes = [axes]

    # Determine global Y-axis limits if not provided
    if y_limits is None:
        max_dev = max(max(values) for positions, values in genome_data.values())
        min_dev = min(min(values) for positions, values in genome_data.values())
    else:
        min_dev, max_dev = y_limits

    for ax, (contig_id, (positions, values)) in zip(axes, genome_data.items()):
        # Plot the data with color-coded bars
        colors = ['#ca0020' if value > 0 else '#0571b0' for value in values]
        ax.bar(positions, values, color=colors, width=step_size / 1e6, label=contig_id)  # Convert step_size to Mb

        # Set contig labels on the left side
        ax.set_ylabel(contig_id, rotation=0, labelpad=20, verticalalignment='center')

        # Set consistent Y-axis limits
        ax.set_ylim(min_dev, max_dev)

        # Move only the ticks to the right side
        ax.yaxis.tick_right()

        # Highlight regions from BED file if provided
        if bed_data is not None and strain_species_id in bed_data and contig_id in bed_data[strain_species_id]:
            for region in bed_data[strain_species_id][contig_id]:
                ax.axvspan(region[0] / 1e6, region[1] / 1e6, color=region[3], alpha=0.5)  # Use specific color code

    # Set X-axis label and figure title
    plt.xlabel("Genomic position (Mb)")
    fig.suptitle(f"{strain_species_id} - GC content deviation")

    # Add Y-axis label on the right side of the figure closer to the plot area
    fig.text(0.95, 0.5, "GC deviation from the mean", va='center', ha='center', rotation='vertical')

    # Save plot as PDF
    plt.savefig(output_file)
    plt.close()

def parse_bed_file(bed_file):
    """Parse BED file and return regions as a dictionary with strain/species ID as a key."""
    bed_data = {}
    with open(bed_file, 'r') as file:
        for line in file:
            columns = line.strip().split()
            strain_species_id = columns[0]
            contig_id = columns[1]
            start = int(columns[2])
            end = int(columns[3])
            region_id = columns[4]
            color_code = columns[5]
            if strain_species_id not in bed_data:
                bed_data[strain_species_id] = {}
            if contig_id not in bed_data[strain_species_id]:
                bed_data[strain_species_id][contig_id] = []
            bed_data[strain_species_id][contig_id].append((start, end, region_id, color_code))
    print(f"Parsed BED data: {bed_data}")
    return bed_data

def process_genome(fasta_file, window_size, step_size, bed_data, output_dir):
    """Process a single genome, calculate GC content and deviation, and generate plots."""
    strain_species_id = os.path.splitext(fasta_file)[0].replace('.scaffolds', '')
    genome_data_gc = {}
    genome_data_gc_deviation = {}
    all_data = []

    with open(fasta_file, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            contig_id = record.id
            sequence = str(record.seq)

            # Calculate mean GC content of the whole genome
            mean_gc_content = calculate_gc_content(sequence)

            # Calculate GC content and positions
            gc_content, positions = sliding_window_gc(sequence, window_size, step_size)

            # Calculate GC content deviation from the mean
            gc_deviation = [gc - mean_gc_content for gc in gc_content]

            # Store data for consolidated plotting
            genome_data_gc[contig_id] = (positions, gc_content)
            genome_data_gc_deviation[contig_id] = (positions, gc_deviation)

            # Append data for consolidated raw data file
            contig_data = pd.DataFrame({'Contig_ID': contig_id, 'Position': positions, 'GC_Content': gc_content, 'GC_Deviation': gc_deviation})
            all_data.append(contig_data)

    # Create output directory for the genome
    genome_output_dir = os.path.join(output_dir, strain_species_id)
    os.makedirs(genome_output_dir, exist_ok=True)

    # Determine global Y-axis limits for non-zoomed plots
    max_gc = max(max(values) for positions, values in genome_data_gc.values())
    min_gc = min(min(values) for positions, values in genome_data_gc.values())
    max_dev = max(max(values) for positions, values in genome_data_gc_deviation.values())
    min_dev = min(min(values) for positions, values in genome_data_gc_deviation.values())

    # Plot consolidated GC content
    plot_consolidated_gc(genome_data_gc, os.path.join(genome_output_dir, f"{strain_species_id}_consolidated_gc_content.pdf"), "GC Content (%)", bed_data, strain_species_id, y_limits=(min_gc, max_gc))

    # Plot consolidated GC content deviation
    plot_consolidated_gc_deviation(genome_data_gc_deviation, os.path.join(genome_output_dir, f"{strain_species_id}_consolidated_gc_deviation.pdf"), step_size, bed_data, strain_species_id, y_limits=(min_dev, max_dev))

    # Save consolidated raw data
    all_data_df = pd.concat(all_data, ignore_index=True)
    raw_table_file = os.path.join(genome_output_dir, f"{strain_species_id}_gc_content.csv")
    all_data_df.to_csv(raw_table_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="GC Content Analysis")
    parser.add_argument('--window_size', type=int, required=True, help="Size of the sliding window")
    parser.add_argument('--step_size', type=int, required=True, help="Step size for sliding window")
    parser.add_argument('--coordinates_file', type=str, help="File with specific genome coordinates for zoomed-in views")
    parser.add_argument('--zoom_window_size', type=int, help="Size of the sliding window for zoomed-in views")
    parser.add_argument('--zoom_step_size', type=int, help="Step size for sliding window for zoomed-in views")
    parser.add_argument('--bed_file', type=str, help="File with regions to highlight")
    parser.add_argument('--output_dir', type=str, default="output", help="Directory for output files")
    parser.add_argument('--num_cpus', type=int, default=os.cpu_count() // 2, help="Number of CPUs to use for parallel processing (default: half of available CPUs)")
    args = parser.parse_args()

    window_size = args.window_size
    step_size = args.step_size
    zoom_window_size = args.zoom_window_size if args.zoom_window_size else window_size
    zoom_step_size = args.zoom_step_size if args.zoom_step_size else step_size
    coordinates_file = args.coordinates_file
    bed_file = args.bed_file
    output_dir = args.output_dir
    num_cpus = args.num_cpus

    bed_data = parse_bed_file(bed_file) if bed_file else None

    # Create the main output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Find all FASTA files in the current directory
    fasta_files = glob.glob("*.fa")  # Assuming FASTA files have .fa extension

    # Process each genome in parallel
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(process_genome, fasta_file, window_size, step_size, bed_data, output_dir)
            for fasta_file in fasta_files
        ]

    # Process zoomed-in views if coordinates file is provided
    if coordinates_file:
        coord_df = pd.read_csv(coordinates_file, sep="\t", header=None, names=["Strain_Species_ID", "Contig_ID", "Start", "End"])
        for index, row in coord_df.iterrows():
            strain_species_id = row["Strain_Species_ID"]
            contig_id = row["Contig_ID"]
            start = int(row["Start"])
            end = int(row["End"])
            fasta_file = f"{strain_species_id}.scaffolds.fa"  # Match the file name format
            with open(fasta_file, 'r') as file:
                for record in SeqIO.parse(file, 'fasta'):
                    if record.id == contig_id:
                        sequence = str(record.seq[start:end])

                        # Calculate mean GC content of the specified region
                        mean_gc_content = calculate_gc_content(sequence)

                        # Calculate GC content and positions
                        gc_content, positions = sliding_window_gc(sequence, zoom_window_size, zoom_step_size, start=start)

                        # Calculate GC content deviation from the mean
                        gc_deviation = [gc - mean_gc_content for gc in gc_content]

                        # Plot GC content
                        zoom_output_dir = os.path.join(output_dir, strain_species_id)
                        os.makedirs(zoom_output_dir, exist_ok=True)
                        plot_consolidated_gc(
                            {contig_id: (positions, gc_content)},
                            os.path.join(zoom_output_dir, f"{strain_species_id}_{contig_id}_{start}_{end}_gc_content.pdf"),
                            "GC Content (%)",
                            bed_data,
                            strain_species_id,
                            figsize=(15, 2)  # Change 2 to your desired vertical size
                        )

                        # Plot GC content deviation
                        plot_consolidated_gc_deviation(
                            {contig_id: (positions, gc_deviation)},
                            os.path.join(zoom_output_dir, f"{strain_species_id}_{contig_id}_{start}_{end}_gc_deviation.pdf"),
                            zoom_step_size,
                            bed_data,
                            strain_species_id,
                            figsize=(15, 2)  # Change 2 to your desired vertical size
                        )

if __name__ == "__main__":
    main()
