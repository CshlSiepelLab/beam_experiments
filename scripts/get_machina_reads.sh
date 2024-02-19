#!/bin/bash

source_directory="/home/staklins/machina/data/sims/m5"
destination_directory="/home/staklins/bayesian_phylogenetic_metastasis/machina_m5_sim_data"

# source_directory="/home/staklins/machina/data/sims/m5"
# destination_directory="/home/staklins/bayesian_phylogenetic_metastasis/machina_m5_sim_data"

# Find all "reads_*.tsv" files in the source directory and its subdirectories
find "$source_directory" -type f -name "reads_*.tsv" -print0 | while IFS= read -r -d $'\0' file; do
    # Extract the wildcard value from the filename
    wildcard_value=$(echo "$file" | grep -oP 'reads_(.*?).tsv' | cut -d '_' -f 2 | cut -d '.' -f 1)

    # Construct the destination directory based on the wildcard value
    destination_subdirectory="$destination_directory/$wildcard_value"

    # # Create the destination directory if it doesn't exist
    # mkdir -p "$destination_subdirectory"

    # Copy the file to the destination directory
    # echo $file
    # echo $destination_subdirectory
    cp "$file" "$destination_subdirectory/"
done