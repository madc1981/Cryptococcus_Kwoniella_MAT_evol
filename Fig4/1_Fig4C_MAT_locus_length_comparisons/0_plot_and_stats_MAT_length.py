#!/usr/bin/env python3

"""
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2024-07-11
Version: 1.0
Description: This script analyzes MAT locus length data and performs statistical comparisons.
             It generates a box plot and a summary table with statistical test results.
Requirements: pandas, matplotlib, seaborn, numpy, scipy, sys, os
Usage: python script_name.py <input_file.tsv>
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu
import sys
import os

def main(input_file):
    # Load the data
    if not os.path.isfile(input_file):
        print(f"Error: The file '{input_file}' does not exist.")
        sys.exit(1)

    data = pd.read_csv(input_file, sep='\t')

    # Convert MAT_length to numeric
    data['MAT_length'] = pd.to_numeric(data['MAT_length'], errors='coerce')

    # Perform the Mann-Whitney U test
    group_alpha = data[data['Mating_Type'] == 'MATalpha']['MAT_length']
    group_a = data[data['Mating_Type'] == 'MATa']['MAT_length']
    stat, p_value = mannwhitneyu(group_alpha, group_a)

    # Calculate means, medians, and standard deviations
    mean_alpha = group_alpha.mean()
    mean_a = group_a.mean()
    median_alpha = group_alpha.median()
    median_a = group_a.median()
    std_alpha = group_alpha.std()
    std_a = group_a.std()

    # Set up plot parameters
    plt.rcParams['pdf.fonttype'] = 42
    plt.figure(figsize=(4, 6))
    palette = {'MATa': '#3cb1d0', 'MATalpha': '#b5d9e8'}
    sns.boxplot(x='Mating_Type', y='MAT_length', data=data, palette=palette, showmeans=True, meanline=True, width=0.4,
                meanprops={"color": "red", "ls": "--", "linewidth": 2},
                medianprops={"color": "blue", "linewidth": 2})

    # Add jittered stripplot to show all data points
    sns.stripplot(x='Mating_Type', y='MAT_length', data=data, color='white', edgecolor='black', linewidth=1, jitter=True, size=8, alpha=0.6)

    # Check if values are larger than 100,000 bp and adjust y-axis labels
    if data['MAT_length'].max() > 100000:
        plt.ylabel('Length (kb)', size=14)
        # Convert values to kb for display
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x/1000:.0f}'))
    else:
        plt.ylabel('Length (bp)', size=14)

    # Add title and labels
    plt.xlabel('Mating Type', size=14)

    # Set font size for x-axis labels
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Add statistical test result at the top with a connecting bar
    significance = 'n.s.' if p_value >= 0.05 else f'p = {p_value:.4g}'
    y_max = data['MAT_length'].max()
    y_range = y_max - data['MAT_length'].min()
    y_bar_height = y_max + 0.05 * y_range
    y_text_height = y_max + 0.052 * y_range
    plt.hlines(y=y_bar_height, xmin=0, xmax=1, color='black', linewidth=1.5)
    plt.text(0.5, y_text_height, significance, ha='center', va='bottom', color='black', fontsize=12)

    # Adjust layout to make room for text
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Generate output file names based on input file name
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    pdf_output_file = f'{base_name}_length_comparison_plot.pdf'
    output_file_path = f'{base_name}_length_comparison_table.tsv'

    # Save the plot as a PDF
    plt.savefig(pdf_output_file, format='pdf')

    # Show the plot
    plt.show()

    # Export the data and statistical test results as a table
    data[['Species_strain_ID', 'Mating_Type', 'MAT_length']].to_csv(output_file_path, sep='\t', index=False)

    with open(output_file_path, 'a') as f:
        f.write(f'\nMann-Whitney U Test Statistic\t{stat}\n')
        f.write(f'Mann-Whitney U Test P-value\t{p_value:.4g}\n')
        f.write(f'\nMean MAT_length (MATalpha)\t{mean_alpha:.1f}\n')
        f.write(f'Standard Deviation MAT_length (MATalpha)\t{std_alpha:.1f}\n')
        f.write(f'Median MAT_length (MATalpha)\t{median_alpha:.1f}\n')
        f.write(f'\nMean MAT_length (MATa)\t{mean_a:.1f}\n')
        f.write(f'Standard Deviation MAT_length (MATa)\t{std_a:.1f}\n')
        f.write(f'Median MAT_length (MATa)\t{median_a:.1f}\n')

    # Print the statistical test result
    print(f'Mann-Whitney U Test:\nStatistic: {stat}\nP-value: {p_value}')
    print(f"Data exported to: {output_file_path}")
    print(f"Plot saved as: {pdf_output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <input_file.tsv>")
        sys.exit(1)
    main(sys.argv[1])

