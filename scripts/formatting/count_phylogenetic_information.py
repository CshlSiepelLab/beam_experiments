
import os
import sys
import pandas as pd
import numpy as np

def count_informative_characters(site_values):
    """
    Count phylogenetically informative characters in a site.
    A character is informative if:
    1. It's not 0 (unedited) or -1 (missing data)
    2. It appears in more than one cell and less than all cells (provides phylogenetic information)
    Returns the number of unique characters that appear more than once.
    """
    # Remove 0s and -1s
    informative_values = site_values[~site_values.isin([0, -1])]
    
    # Count unique values and their frequencies
    value_counts = informative_values.value_counts()
    
    # Count how many unique values appear more than once
    informative_count = len(value_counts[(value_counts > 1) & (value_counts < len(site_values))])
    
    return informative_count



input_file = sys.argv[1]


ext = os.path.splitext(input_file)[1].lower()

if ext == ".tsv":
    sep = "\t"
else:
    sep = ","
    
df = pd.read_csv(input_file, sep=sep, index_col=0)

# Calculate informative characters per site
informative_counts = df.apply(count_informative_characters)

# Calculate ratio of informative characters to total cells
num_informative = sum(informative_counts)
num_cells = df.shape[0]
if num_informative == 0:
    avg_informative = 0.0
else:
    avg_informative = num_informative / num_cells

print(f"Result: {avg_informative:.2f}")
