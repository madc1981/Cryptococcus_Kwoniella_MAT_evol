#!/usr/bin/env python3
"""
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2024-12-20
Version: 2.4
Description:
    This script calculates the GC content of a specified genomic region and compares it to random regions
    either on the same chromosome or across the whole genome. Includes a permutation test for robust
    significance testing. Users can specify the number of random samples and permutations.
Requirements:
    - Biopython
    - SciPy
    - NumPy
Usage:
    python gc_content_analysis.py genome_file region --comparison_mode vs_chr|vs_wg [--excluded_regions excluded_regions.txt] 
        [--num_samples 100] [--num_permutations 10000] output_file.tsv raw_data.tsv
"""

import sys
import argparse
import numpy as np
from Bio import SeqIO
from scipy.stats import ttest_1samp
from random import randint, seed

def parse_excluded_regions(file_path):
    """Parse excluded regions from a tab-separated file."""
    excluded = {}
    if not file_path:
        return excluded
    with open(file_path) as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            if chrom not in excluded:
                excluded[chrom] = []
            excluded[chrom].append((int(start), int(end)))
    return excluded

def is_in_excluded_regions(start, length, excluded_regions, chromosome):
    """Check if a region overlaps any excluded region."""
    end = start + length
    if chromosome not in excluded_regions:
        return False
    for ex_start, ex_end in excluded_regions[chromosome]:
        if start < ex_end and end > ex_start:
            return True
    return False

def random_region_gc(sequence, region_length, excluded_regions, chromosome=None, num_samples=100):
    """Sample GC content from random regions and return their locations and values."""
    seed(42)  # Fixed seed for reproducibility
    sampled_regions = []
    for _ in range(num_samples):
        while True:
            start = randint(0, len(sequence) - region_length)
            if chromosome and is_in_excluded_regions(start, region_length, excluded_regions, chromosome):
                continue
            break
        sample_seq = sequence[start:start + region_length]
        gc_content = calculate_gc(sample_seq)
        sampled_regions.append((start, start + region_length, gc_content))
    return sampled_regions

def calculate_gc(seq):
    """Calculate GC content of a sequence."""
    if len(seq) == 0:
        raise ValueError("Sequence length is zero; unable to calculate GC content.")
    gc_count = sum(1 for base in seq if base in "GCgc")
    return gc_count / len(seq) * 100

def permutation_test(region_gc, sampled_gc, num_permutations=10000):
    """
    Perform a permutation test to assess the significance of GC content differences.
    """
    combined = np.array(sampled_gc + [region_gc])
    observed_diff = abs(region_gc - np.mean(sampled_gc))
    count_extreme = 0

    for _ in range(num_permutations):
        np.random.shuffle(combined)
        perm_region_gc = combined[-1]  # Last value as "region"
        perm_sampled_gc = combined[:-1]  # Rest as "sampled"
        perm_diff = abs(perm_region_gc - np.mean(perm_sampled_gc))
        if perm_diff >= observed_diff:
            count_extreme += 1

    return count_extreme / num_permutations

def main():
    parser = argparse.ArgumentParser(description="GC Content Analysis Script")
    parser.add_argument("genome_file", help="Path to the genome file in GenBank format")
    parser.add_argument("region", help="Region of interest (e.g., chr_2:2955256-2986592)")
    parser.add_argument("--comparison_mode", choices=["vs_chr", "vs_wg"], required=True,
                        help="Comparison mode: vs_chr (same chromosome) or vs_wg (whole genome)")
    parser.add_argument("--excluded_regions", help="Optional file with excluded regions (tab-separated)")
    parser.add_argument("--num_samples", type=int, default=100, help="Number of random regions to sample (default: 100)")
    parser.add_argument("--num_permutations", type=int, default=10000, help="Number of permutations for permutation test (default: 10000)")
    parser.add_argument("output_file", help="File to save statistical results")
    parser.add_argument("raw_data_file", help="File to save raw sampled GC content values")
    args = parser.parse_args()

    # Parse region
    contig_name, coords = args.region.split(":")
    start, end = map(int, coords.split("-"))
    region_length = end - start

    # Parse genome and find the contig of interest
    genome_sequence = ""
    chromosome_sequence = None
    for record in SeqIO.parse(args.genome_file, "genbank"):
        if record.id == contig_name:
            chromosome_sequence = str(record.seq)
        genome_sequence += str(record.seq)
    
    if not chromosome_sequence:
        sys.exit(f"Contig {contig_name} not found in genome file.")

    # Parse excluded regions
    excluded_regions = parse_excluded_regions(args.excluded_regions)

    # Calculate GC content of the specified region and genome
    region_gc_content = calculate_gc(chromosome_sequence[start:end])
    genome_gc_content = calculate_gc(genome_sequence)

    # Sample random regions GC content based on comparison mode
    if args.comparison_mode == "vs_chr":
        sampled_gc = random_region_gc(
            chromosome_sequence,
            region_length,
            excluded_regions,
            chromosome=contig_name,
            num_samples=args.num_samples
        )
    elif args.comparison_mode == "vs_wg":
        sampled_gc = random_region_gc(
            genome_sequence,
            region_length,
            excluded_regions,
            num_samples=args.num_samples
        )

    # Perform t-test and calculate Z-score
    t_stat, ttest_p_value = ttest_1samp([gc for _, _, gc in sampled_gc], region_gc_content)
    sample_mean = np.mean([gc for _, _, gc in sampled_gc])
    sample_std = np.std([gc for _, _, gc in sampled_gc])
    z_score = (region_gc_content - sample_mean) / sample_std

    # Perform permutation test
    perm_p_value = permutation_test(region_gc_content, [gc for _, _, gc in sampled_gc], num_permutations=args.num_permutations)

    # Save raw data to a separate file
    with open(args.raw_data_file, "w") as raw_out:
        raw_out.write("Sample_Start\tSample_End\tGC_Content(%)\n")
        for start, end, gc in sampled_gc:
            raw_out.write(f"{start}\t{end}\t{gc:.2f}\n")

    # Save results to TSV file
    with open(args.output_file, "w") as out:
        out.write("Metric\tValue\n")
        out.write(f"Genome average GC content (%)\t{genome_gc_content:.2f}\n")
        out.write(f"Specified region GC content (%)\t{region_gc_content:.2f}\n")
        out.write(f"Mean GC content of sampled regions (%)\t{sample_mean:.2f}\n")
        out.write(f"T-test statistic\t{t_stat:.2f}\n")
        out.write(f"T-test p-value\t{ttest_p_value:.5f}\n")
        out.write(f"Z-score of specified region\t{z_score:.2f}\n")
        out.write(f"Permutation test p-value\t{perm_p_value:.5f}\n")

    # Print results to console
    print(f"Genome average GC content: {genome_gc_content:.2f}%")
    print(f"Specified region GC content: {region_gc_content:.2f}%")
    print(f"Mean sampled regions GC content: {sample_mean:.2f}%")
    print(f"T-statistic: {t_stat:.2f}, T-test p-value: {ttest_p_value:.5f}")
    print(f"Z-score of specified region: {z_score:.2f}")
    print(f"Permutation test p-value: {perm_p_value:.5f}")
    print(f"Raw data saved to {args.raw_data_file}")

if __name__ == "__main__":
    main()
