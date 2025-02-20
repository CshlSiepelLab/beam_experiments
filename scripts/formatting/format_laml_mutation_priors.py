#!/usr/bin/env python3

import sys

sim_matrix = sys.argv[1]
mutation_priors = sys.argv[2]
outfile = sys.argv[3]

priors = {}
with open(mutation_priors, "r") as f:
    header = f.readline() # skip the header line
    for line in f.readlines():
        code, rate = line.strip().split(",")
        priors[code] = rate


laml_priors = {}
with open(sim_matrix, "r") as f:
    site_numbers = f.readline().strip().split(",")[1:] # the first value is empty at the top left corner of the input sim matrix, above the index cell names
    for line in f.readlines():
        site_codes = line.strip().split(",")[1:] # the first value in the matrix is the cell name, not a mutation
        for site_num, code in enumerate(site_codes):
            if not code.isdigit() or int(code) <= 0:    # skip missing data as "-" or "?" or "-1" and the unedited 0 state
                continue
            if site_num not in laml_priors:
                laml_priors[site_num] = {}
            if code not in laml_priors[site_num]:
                try:
                    laml_priors[site_num][code] = priors[code]
                except KeyError:
                    sys.stderr.write(f"Error: Code {code} not found in priors.\n")
                    sys.exit(1)

with open(outfile, "w") as f:
    for site, code_dict in laml_priors.items():
        for code, rate in code_dict.items():
            f.write(f"{site},{code},{rate}\n")



