#! /usr/bin/env python3

import sys
import numpy as np
import pandas as pd

def minimum_migration_inducing_triplet_count(char_matrix, tissue_labels, primary_tissue):
    """
    Minimum Migration-Inducing Triplet Count (MMITC) with primary tissue exclusion.
    
    Counts pairs of cells (A, B) from different tissues (excluding primary tissue) such that:
      - They share at least one mutation character,
      - The shared mutation character is NOT found in any cell from the primary tissue,
      - They differ by at least one mutation character.
    """
    num_cells = char_matrix.shape[0]
    count = 0
    
    tissues = tissue_labels.values
    cells = char_matrix.index.to_list()
    char_array = char_matrix.to_numpy()
    
    # Extract all mutation states for primary tissue cells
    primary_cells_mask = tissue_labels == primary_tissue
    primary_mutations = char_array[primary_cells_mask.values, :]
    
    # Create a set of mutation states present in primary tissue per site
    # This will be used to exclude shared mutations present in primary tissue
    primary_mutation_sets = [set(primary_mutations[:, col]) for col in range(char_array.shape[1])]
    
    for i in range(num_cells):
        for j in range(i + 1, num_cells):
            # Skip if either cell is from primary tissue or if tissues are same
            if tissues[i] == primary_tissue or tissues[j] == primary_tissue:
                continue
            if tissues[i] == tissues[j]:
                continue
            
            # Positions where mutations are shared
            shared_positions = [k for k in range(char_array.shape[1]) if char_array[i,k] == char_array[j,k]]
            if not shared_positions:
                continue
            
            # Check if at least one shared mutation NOT found in primary tissue at that position
            shared_exclusive = False
            for pos in shared_positions:
                shared_state = char_array[i, pos]
                # If shared state is not in primary tissue mutation set at this site
                if shared_state not in primary_mutation_sets[pos]:
                    shared_exclusive = True
                    break
            if not shared_exclusive:
                continue
            
            # Check if they differ at least one mutation site
            diff = np.sum(char_array[i] != char_array[j])
            if diff < 1:
                continue
            
            count += 1
                
    return count

if __name__ == "__main__":
    # Get input file paths and primary tissue label from command line arguments
    char_matrix_fp = sys.argv[1]  # Path to character matrix file (mutations x cells)
    tissue_labels_fp = sys.argv[2]  # Path to tissue labels file
    primary_tissue = sys.argv[3]  # Label of the primary tissue

    # Load the character matrix (mutations x cells) from TSV file
    char_matrix = pd.read_csv(char_matrix_fp, sep='\t', index_col=0)

    # Load tissue labels and create a mapping from cell names to tissue labels
    tissue_df = pd.read_csv(tissue_labels_fp, sep='\t')
    tissue_dict = pd.Series(tissue_df.tissues.values, index=tissue_df.group_name).to_dict()
    
    # Check for cells that don't have tissue labels and warn if any are found
    missing_cells = [cell for cell in char_matrix.index if cell not in tissue_dict]
    if missing_cells:
        print("Warning: Some cells missing tissue labels, will be ignored:", missing_cells)
    
    # Filter the character matrix to only include cells with known tissue labels
    filtered_cells = [cell for cell in char_matrix.index if cell in tissue_dict]
    char_matrix = char_matrix.loc[filtered_cells]
    
    # Create a Series of tissue labels aligned with the filtered character matrix
    tissues_aligned = pd.Series([tissue_dict[cell] for cell in filtered_cells], index=filtered_cells)
    
    # Calculate and print the Minimum Migration-Inducing Triplet Count
    result = minimum_migration_inducing_triplet_count(char_matrix, tissues_aligned, primary_tissue)
    print("Minimum Migration-Inducing Triplet Count:", result)
