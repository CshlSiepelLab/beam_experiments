
import sys

sim_matrix = sys.argv[1]
mutation_priors = sys.argv[2]
outfile = sys.argv[3]

priors = {}
with open(mutation_priors, "r") as f:
    # skip the header line
    header = f.readline()
    for line in f.readlines():
        code, rate = line.strip().split(",")
        priors[code] = rate


laml_priors = {}
with open(sim_matrix, "r") as f:
    # the first value is empty at the top left corner of the input sim matrix, above the index cell names
    site_numbers = f.readline().strip().split(",")[1:]
    for line in f.readlines():
        # the first value in the matrix is the cell name, not a mutation
        site_codes = line.strip().split(",")[1:]
        for i, code in enumerate(site_codes):
            site_name = i
            # skip missing data as "-" or "?" or "-1" and the unedited 0 state
            if not code.isdigit() or int(code) <= 0:
                # laml still expects the site to be in the output file, even if the code is 0 for all cells
                if site_name not in laml_priors:
                    laml_priors[site_name] = {}
                continue
            # laml wants just int site names for the input priors
            if site_name not in laml_priors:
                laml_priors[site_name] = {}
            if code not in laml_priors[site_name]:
                try:
                    laml_priors[site_name][code] = priors[code]
                except KeyError:
                    sys.stderr.write(f"Error: Code {code} not found in priors.\n")
                    sys.exit(1)

laml_priors = {
    k: v for k, v in sorted(laml_priors.items(), key=lambda item: int(item[0]))
}

# normalize for each site since this is how laml example priors are formatted
for site, code_dict in laml_priors.items():
    total_rate = sum([float(rate) for rate in code_dict.values()])
    for code, rate in code_dict.items():
        laml_priors[site][code] = float(rate) / total_rate

with open(outfile, "w") as f:
    for site, code_dict in laml_priors.items():
        if code_dict:
            for code, rate in code_dict.items():
                f.write(f"{site},{code},{rate}\n")
        # this is again just to satisfy laml requiring all sites to be in the priors even if the site has no mutations
        # alternatively, these sites can be removed from the mutation matrix since they are uninformative
        else:
            f.write(f"{site},0,1\n")
