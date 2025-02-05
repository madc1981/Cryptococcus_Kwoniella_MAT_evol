#!/usr/bin/env python3
"""
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2024-12-20
Version: 1.4
Description:
    This script analyzes codon usage frequencies in the P/R locus and genome-wide coding sequences
    for specified regions in GenBank files. It outputs raw codon-specific data, Chi2 test results, 
    and AT/GC content comparisons. The output filenames are automatically generated using a user-provided prefix.
Requirements:
    - Python >3.7
    - Biopython >1.78
    - SciPy >1.5.2
    - Matplotlib >3.3.2 
    - Seaborn >0.11.0
Usage:
    python 3_codon_usage_analysis.py region_file genbank_folder output_prefix [--no-plot]
"""

import os
import sys
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency

def parse_region_file(region_file):
    """Parse the region file into a dictionary."""
    regions = {}
    with open(region_file, "r") as f:
        for line in f:
            try:
                strain, region = line.strip().split()
                contig, coords = region.split(":")
                start, end = map(int, coords.split("-"))
                regions[strain] = (contig, start, end)
            except ValueError:
                print(f"Warning: Skipping invalid line in region file: {line.strip()}")
    return regions

def extract_cds(genbank_file, contig, start, end):
    """Extract CDS codons within the specified region."""
    pr_cds_codons = []
    genome_cds_codons = []

    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                cds_seq = feature.location.extract(record.seq)
                if len(cds_seq) % 3 != 0:
                    cds_seq = cds_seq[:len(cds_seq) - len(cds_seq) % 3]
                codons = [str(cds_seq[i:i + 3]) for i in range(0, len(cds_seq), 3)]
                genome_cds_codons.extend(codons)
                seq_start = feature.location.start.position
                seq_end = feature.location.end.position
                if record.id == contig and seq_start < end and seq_end > start:
                    pr_cds_codons.extend(codons)

    return pr_cds_codons, genome_cds_codons

def calculate_codon_frequencies(codons):
    """Calculate codon usage frequencies."""
    codon_counter = Counter(codons)
    total_codons = sum(codon_counter.values())
    codon_frequencies = {codon: count / total_codons for codon, count in codon_counter.items()}
    return codon_frequencies, codon_counter

def calculate_at_gc_counts(codons):
    """Calculate AT and GC codon counts."""
    at_count = sum(1 for codon in codons if any(base in codon for base in ["A", "T"]))
    gc_count = sum(1 for codon in codons if any(base in codon for base in ["G", "C"]))
    return at_count, gc_count

def chi_square_test(pr_counter, genome_counter):
    """Perform Chi-Square Test for codon usage."""
    all_codons = set(pr_counter.keys()).union(genome_counter.keys())
    pr_counts = [pr_counter.get(codon, 0) for codon in all_codons]
    genome_counts = [genome_counter.get(codon, 0) for codon in all_codons]
    chi2, p, _, _ = chi2_contingency([pr_counts, genome_counts])
    return chi2, p

def chi_square_at_gc(pr_at_gc, genome_at_gc):
    """Perform a Chi-Square Test on AT/GC counts."""
    chi2, p_value, _, _ = chi2_contingency([pr_at_gc, genome_at_gc])
    return chi2, p_value

def plot_at_gc_comparison(pr_counts, genome_counts, strain):
    """Plot AT/GC Proportions for P/R and Genome-Wide Regions."""
    labels = ["AT", "GC"]
    pr_proportions = [pr_counts[0] / sum(pr_counts), pr_counts[1] / sum(pr_counts)]
    genome_proportions = [genome_counts[0] / sum(genome_counts), genome_counts[1] / sum(genome_counts)]

    x = range(len(labels))
    plt.figure(figsize=(8, 5))
    plt.bar(x, pr_proportions, width=0.4, label="P/R Region", align="center")
    plt.bar([p + 0.4 for p in x], genome_proportions, width=0.4, label="Genome-Wide", align="center")
    plt.xticks([p + 0.2 for p in x], labels)
    plt.ylabel("Proportion")
    plt.title(f"AT/GC Content Comparison - {strain}")
    plt.legend()
    plt.tight_layout()
    plt.show()

def main():
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python 3_codon_usage_analysis.py region_file genbank_folder output_prefix [--no-plot]")
        sys.exit(1)

    region_file = sys.argv[1]
    genbank_folder = sys.argv[2]
    output_prefix = sys.argv[3]
    plot_enabled = "--no-plot" not in sys.argv

    if not os.path.exists(region_file):
        print(f"Error: Region file '{region_file}' not found.")
        sys.exit(1)

    # Auto-generated output filenames
    raw_output_file = f"{output_prefix}_raw_codon_data.tsv"
    stats_output_file = f"{output_prefix}_chi2_results.tsv"
    at_gc_output_file = f"{output_prefix}_at_gc_content.tsv"

    regions = parse_region_file(region_file)
    raw_results = []
    stats_results = []
    at_gc_results = []

    for strain, (contig, start, end) in regions.items():
        print(f"Processing strain: {strain}")
        genbank_file = os.path.join(genbank_folder, f"{strain}.gbk")
        if not os.path.exists(genbank_file):
            print(f"Error: GenBank file for {strain} not found in {genbank_folder}")
            continue

        pr_cds, genome_cds = extract_cds(genbank_file, contig, start, end)
        pr_freqs, pr_counter = calculate_codon_frequencies(pr_cds)
        genome_freqs, genome_counter = calculate_codon_frequencies(genome_cds)
        chi2, p_value = chi_square_test(pr_counter, genome_counter)
        raw_results.extend([[strain, codon, pr_freqs.get(codon, 0), genome_freqs.get(codon, 0)] for codon in pr_freqs.keys()])
        stats_results.append([strain, chi2, p_value])

        # AT/GC Analysis
        pr_at_gc = calculate_at_gc_counts(pr_cds)
        genome_at_gc = calculate_at_gc_counts(genome_cds)
        at_gc_chi2, at_gc_p = chi_square_at_gc(pr_at_gc, genome_at_gc)
        at_gc_results.append([strain, *pr_at_gc, *genome_at_gc, at_gc_chi2, at_gc_p])

        # Visualizations
        if plot_enabled:
            plot_at_gc_comparison(pr_at_gc, genome_at_gc, strain)

    # Save Results
    with open(raw_output_file, "w") as raw_out:
        raw_out.write("Strain\tCodon\tP/R Frequency\tGenome-Wide Frequency\n")
        for row in raw_results:
            raw_out.write("\t".join(map(str, row)) + "\n")

    with open(stats_output_file, "w") as stats_out:
        stats_out.write("Strain\tChi2\tP-Value\n")
        for row in stats_results:
            stats_out.write("\t".join(map(str, row)) + "\n")

    with open(at_gc_output_file, "w") as at_gc_out:
        at_gc_out.write("Strain\tP/R AT Count\tP/R GC Count\tGenome-Wide AT Count\tGenome-Wide GC Count\tChi2\tP-Value\n")
        for row in at_gc_results:
            at_gc_out.write("\t".join(map(str, row)) + "\n")

    print(f"Raw codon data saved to {raw_output_file}")
    print(f"Chi2 test results saved to {stats_output_file}")
    print(f"AT/GC content analysis saved to {at_gc_output_file}")


if __name__ == "__main__":
    main()
