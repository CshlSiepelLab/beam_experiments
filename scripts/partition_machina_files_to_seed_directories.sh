#!/bin/bash

# Iterate over matched pairs
for file in T_seed*.tree; do
    # Extract the seed name from the filename
    seed=$(echo "$file" | awk -F'_' '{print $2}' | awk -F'.' '{print $1}')
    
    # Create a directory with the seed name if it doesn't exist
    mkdir -p "$seed"
    
    # Move the files into the corresponding directory
    mv "$file" "$seed/"
    mv "${file%.tree}.vertex.labeling" "$seed/"
done
