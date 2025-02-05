#!/bin/bash

# Region file specifying genome base names and their respective regions
REGION_FILE="./clade_D_species_PR_regions.txt"

# Directory containing GenBank files (current directory)
GENBANK_DIR="."

# Output summary file
SUMMARY_FILE="gc_content_summary.tsv"

# Write the header line to the summary file
echo -e "Strain\tGenome average GC content (%)\tSpecified region GC content (%)\tT-test statistic\tT-test p-value\tZ-score of specified region\tPermutation test p-value" > "$SUMMARY_FILE"

# Loop through each line in the region file
while read -r genome_base region; do
    # Find the GenBank file that matches the genome base name
    gbk_file=$(find "$GENBANK_DIR" -type f -name "${genome_base}.gbk")
    
    # Check if the GenBank file exists
    if [ -z "$gbk_file" ]; then
        echo "Error: GenBank file for $genome_base not found in $GENBANK_DIR"
        continue
    fi

    # Define output file names based on genome base name and region
    output_file="${genome_base}_stats.tsv"
    raw_data_file="${genome_base}_raw_data.tsv"

    # Run the script for the current GenBank file and region
    python 1b_gc_content_analysis.py "$gbk_file" "$region" --comparison_mode vs_wg \
        --num_samples 1000 --num_permutations 10000 "$output_file" "$raw_data_file"

    # Extract results from the output file and append them to the summary
    if [ -f "$output_file" ]; then
        strain=$(basename "$gbk_file" .gbk)
        genome_gc=$(grep "Genome average GC content (%)" "$output_file" | cut -f2)
        region_gc=$(grep "Specified region GC content (%)" "$output_file" | cut -f2)
        t_stat=$(grep "T-test statistic" "$output_file" | cut -f2)
        t_pval=$(grep "T-test p-value" "$output_file" | cut -f2)
        z_score=$(grep "Z-score of specified region" "$output_file" | cut -f2)
        perm_pval=$(grep "Permutation test p-value" "$output_file" | cut -f2)

        # If any value is missing, set it to "NA" to avoid empty fields
        genome_gc=${genome_gc:-"NA"}
        region_gc=${region_gc:-"NA"}
        t_stat=${t_stat:-"NA"}
        t_pval=${t_pval:-"NA"}
        z_score=${z_score:-"NA"}
        perm_pval=${perm_pval:-"NA"}

        # Append the row to the summary file
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$strain" "$genome_gc" "$region_gc" "$t_stat" "$t_pval" "$z_score" "$perm_pval" >> "$SUMMARY_FILE"
    fi

    echo "Finished processing $genome_base with region $region"
done < "$REGION_FILE"

echo "Summary compiled into $SUMMARY_FILE"
