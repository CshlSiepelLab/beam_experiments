#!/usr/bin/env python3

import sys

original_rates_file = sys.argv[1]
reordering_dict_file = sys.argv[2]
outfile = sys.argv[3]

# get the original code to rate mapping
original_rates = {}
with open(original_rates_file) as f:
    # skip the header
    f.readline()
    for line in f:
        code, rate = line.strip().split(",")
        original_rates[code] = rate

# get the mapping of new code to old code
new_rates = {}
with open(reordering_dict_file) as f:
    # skip the header
    f.readline()
    for line in f:
        new_code, old_code = line.strip().split(",")
        # we dont need to map missing data code to a rate since beam only takes in fixed rates for mutation outcomes, not silenced sites
        if str(new_code) != "-1":
            new_rates[new_code] = original_rates[old_code]

# order the new_rates by the new codes in ascending order
ordered_rates = [rate for code, rate in sorted(new_rates.items(), key=lambda item: int(item[0]))]

with open(outfile, "w") as f:
    for rate in ordered_rates:
        f.write(f"{rate}\n")
