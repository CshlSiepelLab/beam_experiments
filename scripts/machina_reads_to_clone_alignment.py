#!/usr/bin/env python3

### This script takes in the MACHINA sim data tsv format and makes clone alignment to then use in PathFinder setup to get trees with branch lengths resolved.

import sys
import pandas as pd

reads_tsv=sys.argv[1]

# reads_tsv="machina_data/sims/machina_m5_sim_data/seed0/reads_seed0.tsv"

colnames = ["sample_index", "sample_label", "anatomical_site_index", "anatomical_site_label", "character_index", "character_label", "ref", "var"]
reads_df = pd.read_csv(reads_tsv, sep='\t', names=colnames, skiprows=4)

cols = ["sample_index", "var"]
reads_df_subset = reads_df.loc[:, cols]
grouped_df = reads_df_subset.groupby("sample_index")

alignment = ""
for index, group in grouped_df:
    alignment += f">Clo{index}\n"
    for row in group.iterrows():
        if int(row[1]["var"]) == 0:
            alignment += "T"
        else:
            alignment += "A"
    alignment += "\n"

outfile = reads_tsv.split(".")[0] + ".fas"
with open(outfile, "w") as file:
    file.write(alignment)