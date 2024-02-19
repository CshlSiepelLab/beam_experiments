#!/usr/bin/env python3

### This script takes in the MACHINA sim data tsv format and makes clone alignment to then use in PathFinder setup to get trees with branch lengths resolved.

import sys
import pandas as pd

# reads_tsv=sys.argv[1]

reads_tsv="machina_m5_sim_data/seed0/reads_seed0.tsv"

colnames = ["sample_index", "sample_label", "anatomical_site_index", "anatomical_site_label", "character_index", "character_label", "ref", "var"]
reads_df = pd.read_csv(reads_tsv, sep='\t', names=colnames, skiprows=4)

reads_df