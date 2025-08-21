
import sys
import pandas as pd
import numpy as np

def count_informative_characters(site_values):
    """
    Count phylogenetically informative characters in a site.
    A character is informative if:
    1. It's not 0 (unedited) or -1 (missing data)
    2. It appears in more than one cell (provides phylogenetic information)
    Returns the number of unique characters that appear more than once.
    """
    # Remove 0s and -1s
    informative_values = site_values[~site_values.isin([0, -1])]
    
    # Count unique values and their frequencies
    value_counts = informative_values.value_counts()
    
    # Count how many unique values appear more than once
    informative_count = len(value_counts[value_counts > 1])
    
    return informative_count



input_file = sys.argv[1]

# Read the TSV file
df = pd.read_csv(input_file, sep='\t', index_col=0)

# Calculate informative characters per site
informative_counts = df.apply(count_informative_characters)

# Calculate average informative characters
avg_informative = informative_counts.mean()

# Print results
print("\nPhylogenetically informative characters per site:")
print("-" * 50)
for site, count in informative_counts.items():
    print(f"Site {site}: {count} informative characters")

print("\nSummary Statistics:")
print("-" * 50)
print(f"Average informative characters per site: {avg_informative:.2f}")
print(f"Total number of sites: {len(informative_counts)}")
print(f"Total number of informative characters: {informative_counts.sum()}")