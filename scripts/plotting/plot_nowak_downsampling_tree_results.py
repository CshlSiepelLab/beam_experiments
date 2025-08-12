
import pandas as pd
import sys
import matplotlib.pyplot as plt

csv_path = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24/downsampling_counts.csv"

# Read the CSV file
df = pd.read_csv(csv_path)

# Plot the data
fs = 14

plt.figure(figsize=(12, 6))

# Group the data by 'mouse' and 'cp'
grouped = df.groupby(["mouse", "cp"]).sum().reset_index()

# Create the bar plot
bar_width = 0.35
index = range(len(grouped))

fig, ax = plt.subplots()

bar1 = ax.bar(index, grouped["original_num_asvs"], bar_width, label="Original")
bar2 = ax.bar(
    [i + bar_width for i in index],
    grouped["downsample_num_asvs"],
    bar_width,
    label="Downsampled",
)

# Add labels and title
ax.set_ylabel("Number of tips", fontsize=fs)
ax.set_xticks([i + bar_width / 2 for i in index])
ax.tick_params(axis="x", labelsize=fs)
ax.tick_params(axis="y", labelsize=fs)
ax.set_xticklabels(
    [f"{row['mouse']}-{row['cp']}" for _, row in grouped.iterrows()], rotation=90
)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False, fontsize=fs)

plt.tight_layout()
plt.savefig(csv_path.replace(".csv", ".pdf"))
plt.close()
