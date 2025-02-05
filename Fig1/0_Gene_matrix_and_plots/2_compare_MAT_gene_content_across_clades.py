#!/usr/bin/env python3

"""
Author: Marco A. Coelho @ Heitman Lab, Duke University
Date: 2024-12-01
Version: 1.2
Description: This script compares genes between specified clades based on gene presence/absence data.
             Features:
             - Compares gene presence across multiple clades
             - Option to exclude HD locus-associated genes
             - Generates comprehensive comparison table with totals
             - Creates Venn diagrams for 2-3 clades showing gene set overlaps
Requirements: matplotlib_venn, pandas, matplotlib

Usage:
    python 2_compare_MAT_gene_content_across_clades.py input_file.tsv CladeA CladeB ... [--exclude-hd-genes]
Example:
    python 2_compare_MAT_gene_content_across_clades.py MAT_loci_genes_presence_absence_matrix_with_clades.tsv A B C2 D F G H I
    python 2_compare_MAT_gene_content_across_clades.py MAT_loci_genes_presence_absence_matrix_with_clades.tsv A B C --exclude-hd-genes
"""

import pandas as pd
from itertools import combinations
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
import os
import sys

# Define HD locus genes (adjust as needed)
HD_GENES = ['SXI1', 'SXI2', 'RPL22', 'CAP1', 'SPO14']

def compare_clades_with_totals(input_file, clades, exclude_hd_genes=False):
    """
    Compares gene presence/absence among specified clades, generates a results table with totals, and optionally
    creates Venn diagrams for clade comparisons. Optionally excludes HD-associated genes.
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

    # Load the input dataset
    data = pd.read_csv(input_file, sep='\t')

    # Ensure the clades exist in the dataset
    available_clades = data['Clade'].unique().tolist()
    for clade in clades:
        if clade not in available_clades:
            raise ValueError(f"Clade '{clade}' is not present in the input file.")

    # Include all gene columns for analysis, optionally excluding HD-associated genes
    gene_columns = data.columns[3:]  # Assumes genes start from the 4th column onward
    if exclude_hd_genes:
        gene_columns = [col for col in gene_columns if col not in HD_GENES]

    # Filter and process data for each clade
    clade_gene_sets = {}
    for clade in clades:
        clade_data = data[data['Clade'] == clade].copy()
        clade_data[gene_columns] = clade_data[gene_columns].replace({'pseudo': 1, '?': 0}).astype(int)
        clade_gene_sets[clade] = set(clade_data[gene_columns].columns[clade_data[gene_columns].any(axis=0)])

    # Define Cryptococcus and Kwoniella clades for group-level comparisons
    cryptococcus_clades = ['A', 'B', 'C1', 'C2', 'D']
    kwoniella_clades = ['E', 'F', 'G', 'H', 'I']

    cryptococcus_genes = set.intersection(
        *[clade_gene_sets[clade] for clade in cryptococcus_clades if clade in clade_gene_sets]
    )
    kwoniella_genes = set.intersection(
        *[clade_gene_sets[clade] for clade in kwoniella_clades if clade in clade_gene_sets]
    )
    shared_genes = cryptococcus_genes & kwoniella_genes

    # Compile results into a single table
    all_genes = [gene for gene in gene_columns]
    results = []
    pairwise_combos = list(combinations(clades, 2))
    for gene in all_genes:
        gene_presence = {clade: "Yes" if gene in clade_gene_sets.get(clade, set()) else "No" for clade in clades}
        is_common = "Yes" if all(gene_presence[clade] == "Yes" for clade in clades) else "No"
        is_unique = {clade: "Yes" if gene_presence[clade] == "Yes" and all(
            gene_presence[other_clade] == "No" for other_clade in clades if other_clade != clade) else "No" for clade in clades}
        pairwise_presence = [("Yes" if gene_presence[clade1] == "Yes" and gene_presence[clade2] == "Yes" else "No")
                             for clade1, clade2 in pairwise_combos]
        in_cryptococcus = "Yes" if gene in cryptococcus_genes else "No"
        in_kwoniella = "Yes" if gene in kwoniella_genes else "No"
        common_crypt_kwoniella = "Yes" if gene in shared_genes else "No"
        results.append([gene] + [gene_presence[clade] for clade in clades] + [is_common] + 
                       [is_unique[clade] for clade in clades] + pairwise_presence + 
                       [in_cryptococcus, in_kwoniella, common_crypt_kwoniella])

    output_columns = (
        ["Gene"]
        + [f"In_clade_{clade}" for clade in clades]
        + ["Common"]
        + [f"Unique_to_clade_{clade}" for clade in clades]
        + [f"Pairwise_{clade1}_{clade2}" for clade1, clade2 in pairwise_combos]
        + ["In_Cryptococcus", "In_Kwoniella", "Common_Cryptococcus_Kwoniella"]
    )
    results_df = pd.DataFrame(results, columns=output_columns)

    # Add totals
    totals = results_df.iloc[:, 1:].apply(lambda col: col.value_counts().get("Yes", 0))
    totals_row = pd.DataFrame([[totals.get(col, "Total") if col != "Gene" else "Total" for col in results_df.columns]],
                              columns=results_df.columns)
    results_df_with_totals = pd.concat([results_df, totals_row], ignore_index=True)

    # Modify output filenames based on the --exclude-hd-genes flag
    suffix = "_noHD" if exclude_hd_genes else ""
    output_file = f'comparison_with_totals_{"_".join(clades)}{suffix}.tsv'
    results_df_with_totals.to_csv(output_file, sep='\t', index=False)
    print(f"Results with totals saved to {output_file}")

    # Update Venn diagram filenames as well
    if len(clades) == 2:
        venn_file = f'venn_{"_".join(clades)}{suffix}.pdf'
        plt.figure(figsize=(8, 8))
        venn2([clade_gene_sets[clades[0]], clade_gene_sets[clades[1]]], 
              (f'Clade {clades[0]}', f'Clade {clades[1]}'))
        plt.savefig(venn_file)
        print(f"Venn diagram saved to {venn_file}")
        plt.show()
    elif len(clades) == 3:
        venn_file = f'venn_{"_".join(clades)}{suffix}.pdf'
        plt.figure(figsize=(8, 8))
        venn3([clade_gene_sets[clades[0]], clade_gene_sets[clades[1]], clade_gene_sets[clades[2]]], 
              (f'Clade {clades[0]}', f'Clade {clades[1]}', f'Clade {clades[2]}'))
        plt.savefig(venn_file)
        print(f"Venn diagram saved to {venn_file}")
        plt.show()
    else:
        print("Venn diagrams are only supported for 2 or 3 clades.")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python 2_compare_MAT_gene_content_across_clades.py MAT_loci_genes_presence_absence_matrix_with_clades.tsv Clade1 Clade2 ... [--exclude-hd-genes]")
        sys.exit(1)

    input_file = sys.argv[1]
    clades = [arg for arg in sys.argv[2:] if not arg.startswith("--")]
    exclude_hd_genes = "--exclude-hd-genes" in sys.argv

    compare_clades_with_totals(input_file, clades, exclude_hd_genes)
