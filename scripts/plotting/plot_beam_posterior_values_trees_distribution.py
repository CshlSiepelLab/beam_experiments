#!/usr/bin/env python3

import re
import sys
import os
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import random
import matplotlib.pyplot as plt


posterior_file = "/grid/siepel/home_norepl/staklins/stephen_data/beast_bayesian_migration_graph_inference/snakemake_performance_repeat_origin_scaling_implemented_10_15_24_uniform_50cells_50sites_data_7_24_24/metastabayes/pR_5978/combined.trees"

outfile = posterior_file.replace(".trees", "_posterior_distribution.pdf")

names_dict = {}
trees = []

with open(posterior_file, "r") as file:
    for line in file:
        line = line.strip()
        if line.startswith("tree"):
            trees.append(line)
        # lines that begin with a number and have two fields are translate lines
        elif line and line[0].isdigit() and len(line.split(" ")) > 1:
            key_value = line.split(" ")
            key = key_value[0]
            # remove trailing comma for translate values
            value = key_value[1].replace(",", "")
            names_dict[key] = value

# sort trees by posterior
pattern = re.compile(r"tree STATE_\d+ = \[&posterior=(-?\d+\.\d+),")
sorted_trees = sorted(
    trees, key=lambda s: float(pattern.search(s).group(1)), reverse=True
)
posterior_values = re.findall(r"\[&posterior=(-?\d+\.\d+),", "".join(sorted_trees))
posterior_values = [round(float(value), 2) for value in posterior_values]

# get peak value and density
kde = gaussian_kde(posterior_values)
x_values = np.linspace(min(posterior_values), max(posterior_values), 1000)

# plot posterior values to see peak
fs = 10
plt.plot(x_values, kde(x_values))
plt.hist(posterior_values, bins=100, density=True, alpha=0.5, color="grey")
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel("Posterior value", fontsize=18)
plt.ylabel("Frequency", fontsize=18)
plt.tight_layout()
plt.savefig(outfile)
plt.close()
