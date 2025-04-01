#!/bin/bash

main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k"

# Define the list of files
files=$(find /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k -maxdepth 1 -type f -name "*_all_consensus_classifications.csv")

# Define the output file
output_file="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/all_thresholds_classification_summary.csv"

# Write the header to the output file
echo "consensus_threshold,met_to_met,met_to_primary" > "$output_file"

# Loop through each file
for file in $files; do
    echo "Processing file: $file"
    
    # Extract the threshold from the filename
    threshold=$(basename "$file" | cut -d'_' -f1)
    echo "Threshold: $threshold"
    
    # Calculate the percentage of True values in each column (excluding the first column)
    percentages=$(awk -F, '
    NR > 1 {
        for (i = 2; i <= NF; i++) {
            if ($i == "True") {
                true_count[i]++
            }
            total_count[i]++
        }
    }
    END {
        for (i = 2; i <= NF; i++) {
            percentage = (true_count[i] / total_count[i]) * 100
            printf "%.2f,", percentage
        }
    }
    ' "$file")
    
    # Remove the trailing comma from percentages
    percentages=${percentages%,}
    
    # Write the threshold and percentages to the output file
    echo "$threshold,$percentages" >> "$output_file"
    
    echo
done

