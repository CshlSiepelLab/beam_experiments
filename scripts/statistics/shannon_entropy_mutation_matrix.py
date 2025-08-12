
import sys
import pandas as pd
import scipy.stats
import numpy as np


file_path = sys.argv[1]

# read in mutation matrix
matrix = pd.read_csv(file_path, sep="\t", header=0, index_col=0)

# count the occurence of each mutation across the matrix
probs = {}
total = 0
for row in matrix.iterrows():
    for value in row[1:]:
        values = value.values
        for mut in values:
            if mut not in probs:
                print(mut)
                probs[mut] = 1
            else:
                probs[mut] += 1
            total += 1

total_probs = len(probs)
probs_array = np.zeros(total_probs)

i = 0

for key, value in probs.items():
    probs_array[i] = value / total
    i += 1
scipy_entropy = scipy.stats.entropy(probs_array, base=2)

# print(f"Shannon Entropy manual: {entropy}")
print(f"Shannon Entropy scipy: {scipy_entropy}")
