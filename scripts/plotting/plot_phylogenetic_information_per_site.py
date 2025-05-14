#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


variable_file = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_2_25_25_data_from_8_19_24/phylogenetic_information_per_site_variable_rates.csv"
quinn_file = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/phylogenetic_information_per_site_quinn.csv"
serio_file = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_2_24_25/phylogenetic_information_per_site_serio.csv"
outfile = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_2_24_25/phylogenetic_information_per_site.pdf"


variable_data = pd.read_csv(variable_file)
quinn_data = pd.read_csv(quinn_file)
serio_data = pd.read_csv(serio_file)

# get groupings for variable rate data based on mutation rates
variable_data['Mutation rate'] = variable_data['name'].str.split('_').str[1].str.replace('mut', '0.')
variable_data['Mutation rate'] = variable_data['Mutation rate'].astype(float)


# Create figure with 3 subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
fs=22
num_bins = 100

# Plot for variable
sns.histplot(data=variable_data, 
            x='average_informative_characters_per_site',
            hue='Mutation rate',
            palette='tab10',
            ax=ax1, kde=False, bins=num_bins)
ax1.set_title('Variable rates simulated dataset', fontsize=fs)
ax1.set_xlabel('', fontsize=fs)
ax1.set_ylabel('Count', fontsize=fs)



# Plot for quinn
sns.histplot(data=quinn_data, 
            x='average_informative_characters_per_site',
            ax=ax2, kde=False, bins=num_bins, color='grey')
ax2.set_title('Lung cancer dataset', fontsize=fs)
ax2.set_xlabel('', fontsize=fs)
ax2.set_ylabel('Count', fontsize=fs)

# Plot for serio
sns.histplot(data=serio_data, 
            x='average_informative_characters_per_site',
            ax=ax3, kde=False, bins=num_bins, color='grey')
ax3.set_title('Prostate cancer dataset', fontsize=fs)
ax3.set_xlabel('Mean site count of unique characters\nshared by more than one cell', fontsize=fs)
ax3.set_ylabel('Count', fontsize=fs)

# Get the maximum y value across all plots
max_x = max(ax1.get_xlim()[1], ax2.get_xlim()[1], ax3.get_xlim()[1])

# Set the same y limits for all plots
ax1.set_xlim(0, max_x)
ax2.set_xlim(0, max_x)
ax3.set_xlim(0, max_x)

# ax3.set_ylim(0,10)
# ax3.set_yticks(range(0, 11, 1))

# Remove grid
ax1.grid(False)
ax2.grid(False)
ax3.grid(False)

ax1.tick_params(labelsize=fs)
ax2.tick_params(labelsize=fs)
ax3.tick_params(labelsize=fs)

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()


# in-vitro data
invitro = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/in_vitro_data_4_24_25/phylogenetic_information_per_site_in_vitro_sorted_subset_cp_above_10.csv"
outfile_invitro = invitro.replace(".csv", ".pdf")

invitro_data = pd.read_csv(invitro)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
fs=22
num_bins = 100


# Plot for invitro
sns.histplot(data=invitro_data, 
            x='average_informative_characters_per_site',
            ax=ax1, kde=False, bins=num_bins, color='grey')
ax1.set_title('In-vitro dataset', fontsize=fs)
ax1.set_xlabel('', fontsize=fs)
ax1.set_ylabel('Count', fontsize=fs)
ax1.set_xlabel('Mean site count of unique characters\nshared by more than one cell', fontsize=fs)


# Get the maximum y value across all plots
max_x = 9.358

ax1.set_xlim(0, max_x)
ax1.grid(False)
ax1.tick_params(labelsize=fs)

# in vitro zoom > 1
invitro_data_zoom = invitro_data[invitro_data['average_informative_characters_per_site'] > 1]
sns.histplot(data=invitro_data_zoom, 
            x='average_informative_characters_per_site',
            ax=ax3, kde=False, bins=num_bins, color='grey')
ax3.set_title('In-vitro dataset', fontsize=fs)
ax3.set_xlabel('', fontsize=fs)
ax3.set_ylabel('Count', fontsize=fs)
ax3.set_xlabel('Mean site count of unique characters\nshared by more than one cell', fontsize=fs)
ax3.set_xlim(0, max_x)
ax3.set_ylim(0, 16)
ax3.set_yticks(range(0, 16, 2))

ax3.grid(False)
ax3.tick_params(labelsize=fs)

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig(outfile_invitro, dpi=300, bbox_inches='tight')
plt.close()
