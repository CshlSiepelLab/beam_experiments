#!/usr/bin/env python3

### This script reads in a BEAST posterior file for trees labeld with tissues at each node and then subsets all sampled trees for a specified number of trees with the highest posterior probability and then collapses these to migraiton graphs.

import re, sys
import numpy as np
import networkx as nx

def find_max_bins(data):
    max_bins = 1
    max_bin_freq = 0
    for num_bins in range(2, len(data) + 1):
        hist, _ = np.histogram(data, bins=num_bins)
        max_freq = np.max(hist)
        num_max_freq_bins = np.sum(hist == max_freq)
        if num_max_freq_bins == 1:
            max_bin_freq = max_freq
            max_bins = num_bins
    return max_bins


# posterior_file = sys.argv[1]

posterior_file = "beast_gundem_2015_2_21_24/A10_sym/tissue_tree_with_trait.trees"

# set the number of trees to obtain as graphs
n = 3

names_dict = {}
trees = []

with open(posterior_file, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('tree'):
            trees.append(line)
        # lines that begin with a number are translate lines
        elif line and line[0].isdigit():
            key_value = line.split(' ')
            key = key_value[0]
            # remove trailing comma for translate values
            value = key_value[1][0:-1]
            names_dict[key] = value

# sort trees by posterior
pattern = re.compile(r'tree STATE_\d+ = \[&posterior=(-?\d+\.\d+)\]')
sorted_trees = sorted(trees, key=lambda s: float(pattern.search(s).group(1)), reverse = True)

# get the top n trees with highest posterior
top_n_trees = sorted_trees[0:n]

# get the top n maximum clade credibility trees by finding closest to the peak of probability density function
posterior_values = re.findall(r'\[&posterior=(-?\d+\.\d+)\]', "".join(trees))
posterior_values = [round(float(value), 2) for value in posterior_values]

max_bins = find_max_bins(posterior_values)

