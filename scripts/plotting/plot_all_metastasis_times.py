#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


file_path = sys.argv[1]
origin_time = float(sys.argv[2])
origin_tissue = sys.argv[3]
outfile = sys.argv[4]


# file_path = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/all_expected_migration_times_no_multiedge.csv'
# origin_time = 54
# origin_tissue = "LL"
# outfile = "./test.pdf"


data = pd.read_csv(file_path)

# get the colors for each tissue
DEFAULT_COLORS = ["#006400", "#FF0000", "#0000CD", "#FFA500", "#800080", "#808080", "#FFC0CB", "#ADD8E6", "#A52A2A", "#FFFF00"]*3
all_tissues = list(sorted(set([tis.split("_")[1] for tis in data['source_recipient'].values]) - {origin_tissue}))
custom_colors = {node: color for node, color in zip(all_tissues, DEFAULT_COLORS[0:len(all_tissues)])}
all_tissues = [origin_tissue] + all_tissues
custom_colors[origin_tissue] = "black"

# Calculate the mean values for sorting

data['source_recipient'] = [route.replace("_", "->") for route in data['source_recipient']]

custom_colors_new = {}
for source_recipient in pd.unique(data["source_recipient"]):
    source, recipient = source_recipient.split("->")
    custom_colors_new[source_recipient] = custom_colors[recipient]

# Group by source tissue and then sort by time within each source tissue
data['source'] = data['source_recipient'].apply(lambda x: x.split("->")[0])
data['recipient'] = data['source_recipient'].apply(lambda x: x.split("->")[1])
data['mean_time'] = data.groupby('source_recipient')['mid_time'].transform('mean')
data = data.sort_values(by=['source', 'mean_time'])

# Create the plot
plt.figure(figsize=(12, 6))

sns.boxplot(x='mid_time', y='source_recipient', data=data, order=data['source_recipient'].unique(), orient='h', palette=custom_colors_new, showfliers=False, boxprops=dict(alpha=0.5))
sns.stripplot(x='mid_time', y='source_recipient', data=data, order=data['source_recipient'].unique(), color='black', alpha=0.5, jitter=True)

plt.xlim(0, origin_time)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.xlabel('Time', fontsize=22)
plt.ylabel('Source tissue -> Recipient tissue', fontsize=22)

ax = plt.gca()
for label in ax.get_yticklabels():
    source, recipient = label.get_text().split("->")
    label.set_color(custom_colors[source])
    label.set_text(f"{source}->{recipient}")

plt.tight_layout()
plt.savefig(outfile)
plt.close()


# Create a new plot for migration trajectories
plt.figure(figsize=(12, 6))

# Create a dictionary to map tissues to y-axis positions
tissue_positions = {tissue: i for i, tissue in enumerate(all_tissues)}

# Plot vertical lines for each migration event
for _, row in data.iterrows():
    source = row['source']
    recipient = row['recipient']
    mid_time = row['mid_time']
    plt.arrow(mid_time, tissue_positions[source], 0, tissue_positions[recipient] - tissue_positions[source], 
              color=custom_colors[recipient], linestyle='-', linewidth=2, alpha=0.7, head_width=0.8, head_length=0.125, length_includes_head=True)

# Set y-axis labels and ticks
plt.yticks(range(len(all_tissues)), all_tissues, fontsize=22)
plt.xticks(fontsize=22)
plt.xlabel('Time', fontsize=22)
plt.ylabel('Tissues', fontsize=22)

# Add grid lines for better readability
ax = plt.gca()
for label in ax.get_yticklabels():
    tissue = label.get_text()
    label.set_color(custom_colors[tissue])

# Draw horizontal lines for each tissue
for tissue, position in tissue_positions.items():
    if tissue == origin_tissue:
        plt.axhline(y=position, xmin=0, xmax=origin_time, color=custom_colors[tissue], linestyle='-', linewidth=10, alpha=1.0)
    else:
        first_event_time = data[data['recipient'] == tissue]['mid_time'].min()
        plt.axhline(y=position, xmin=first_event_time/origin_time, xmax=origin_time, color=custom_colors[tissue], linestyle='-', linewidth=10, alpha=1.0)

plt.tight_layout()
plt.savefig(outfile.replace('.pdf', '_trajectories.pdf'))
plt.close()
