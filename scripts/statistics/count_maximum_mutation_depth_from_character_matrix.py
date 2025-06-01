#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from itertools import combinations

def count_shared_mutations(df):
    """
    Compute the number of shared mutations for all pairwise cell comparisons.
    Ignores 0 and -1 values as they represent missing data and unedited states.
    
    Args:
        df: pandas DataFrame with cells as rows and mutation sites as columns
    
    Returns:
        max_shared: maximum number of shared mutations found
        max_pair: tuple of cell names with maximum shared mutations
    """
    cells = df.index.tolist()
    max_shared = 0
    max_pair = None
    
    # Iterate through all possible cell pairs
    for cell1, cell2 in combinations(cells, 2):
        row1 = df.loc[cell1]
        row2 = df.loc[cell2]
        
        # Only consider positions where both cells have valid mutations (not 0 or -1)
        valid_mask = ~((row1.isin([0, -1])) | (row2.isin([0, -1])))
        
        # Count where mutations match and are valid
        shared = np.sum((row1 == row2) & valid_mask)
        
        if shared > max_shared:
            max_shared = shared
            max_pair = (cell1, cell2)
    
    return max_shared, max_pair

def main():
    # Read the TSV file
    input_file = sys.argv[1]

    df = pd.read_csv(input_file, sep='\t', index_col=0)

    # Convert all data to integers
    df = df.astype(int)

    # Compute maximum shared mutations
    max_shared, max_pair = count_shared_mutations(df)
    
    # Print results
    print(f"Maximum shared mutation depth: {max_shared}")

if __name__ == "__main__":
    main() 