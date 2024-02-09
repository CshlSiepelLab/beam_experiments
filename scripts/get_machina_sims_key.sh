#!/bin/bash

main_directory="/local/storage/no-backup/staklins-scratch/machina/data/sims/m8"

# Set the path to your output .tsv file
output_tsv="/local/storage/no-backup/staklins-scratch/bayesian_phylogenetic_metastasis/machina_m8_sim_data/machina_sims_m8_seeding_routes_key.tsv"

# Create the header for the .tsv file
echo -e "Subdirectory\tValues" > "$output_tsv"

# Find all subdirectories in the main directory
find "$main_directory" -type d -mindepth 1 -maxdepth 1 | while read -r subdirectory; do
    # Extract subdirectory name
    subdirectory_name=$(basename "$subdirectory")
    
    for file in ${subdirectory}/T_*.tree;
    do
        file_name=$(basename "$file")
        value=$(echo $file_name | awk -F'_' '{print $2}' | awk -F'.' '{print $1}')
        echo -e "$subdirectory_name\t$value" >> "$output_tsv"
    done
done
