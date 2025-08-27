
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV file
csv_file = "/grid/siepel/home/staklins/stored_results/beam/latest_results/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_24_24/beam_gtr/all_posterior_expected_migration_counts_over_parsimony.csv"

data = pd.read_csv(csv_file)

# Ungroup the data
data = data.reset_index(drop=True)

# Calculate the difference between the two columns
data["difference"] = data["beast_migration_count"] - data["parsimony_migration_count"]

# Plot the distribution of the differences
grouped = data.groupby("name").mean()

# Plot the average difference grouped by sim_name, so compute the mean difference for each group and plot those values
fontsize = 24
plt.figure(figsize=(10, 6))
plt.hist(grouped["difference"], bins=100, edgecolor="black", color="grey")
plt.xlabel(
    "Posterior expected\n[BEAM migration count - Parsimony migration count]",
    fontsize=fontsize,
)
plt.ylabel("Number of datasets", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.xlim(0, 25)
plt.ylim(0, 30)
# plt.xticks(np.arange(0, 5.5, 1), fontsize=fontsize)
# plt.yticks(np.arange(0, 45, 5), fontsize=fontsize)
plt.tight_layout()
outfile = csv_file.replace(".csv", ".pdf")
plt.savefig(outfile)
plt.close()
