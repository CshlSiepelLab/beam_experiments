
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


variable_file = "/grid/siepel/home/staklins/stored_results/beam/latest_results/variable_migration_and_mutation_rates_data_8_19_24/phylogenetic_information_to_cell_ratio_variable_rates.csv"
quinn_file = "/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/phylogenetic_information_to_cell_ratio_quinn.csv"
serio_file = "/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/phylogenetic_information_to_cell_ratio_serio.csv"
simeonov_file = "/grid/siepel/home/staklins/stored_results/beam/latest_results/simeonov_2021_pancreatic_cancer_data/phylogenetic_information_to_cell_ratio_simeonov.csv"
yang_file = "/grid/siepel/home/staklins/stored_results/beam/latest_results/yang_2022_lung_cancer_data/phylogenetic_information_to_cell_ratio_yang.csv"

outfile = "/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/phylogenetic_information_to_cell_ratio.pdf"


variable_data = pd.read_csv(variable_file)
quinn_data = pd.read_csv(quinn_file)
serio_data = pd.read_csv(serio_file)
simeonov_data = pd.read_csv(simeonov_file)
yang_data = pd.read_csv(yang_file)

# get groupings for variable rate data based on mutation rates
variable_data['Mutation rate'] = variable_data['name'].str.split('_').str[1].str.replace('mut', '0.')
variable_data['Mutation rate'] = variable_data['Mutation rate'].astype(float)


# Create figure with subplots
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(10, 20))
fs=22

# Plot for variable
sns.histplot(data=variable_data, 
            x='informative_muts_to_cell_ratio',
            hue='Mutation rate',
            palette='tab10',
            ax=ax1, kde=False, bins=100)
ax1.set_title('Variable rates simulated dataset', fontsize=fs)
ax1.set_xlabel('', fontsize=fs)
ax1.set_ylabel('Count', fontsize=fs)



# Plot for quinn
sns.histplot(data=quinn_data, 
            x='informative_muts_to_cell_ratio',
            ax=ax2, kde=False, bins=30, color='grey')
ax2.set_title('Lung cancer dataset', fontsize=fs)
ax2.set_xlabel('', fontsize=fs)
ax2.set_ylabel('Count', fontsize=fs)

# Plot for serio
sns.histplot(data=serio_data, 
            x='informative_muts_to_cell_ratio',
            ax=ax3, kde=False, bins=25, color='grey')
ax3.set_title('Prostate cancer dataset', fontsize=fs)
ax3.set_xlabel('Ratio of phylogenetically informative mutations to cells', fontsize=fs)
ax3.set_ylabel('Count', fontsize=fs)

# Plot for simeonov
sns.histplot(data=simeonov_data, 
            x='informative_muts_to_cell_ratio',
            ax=ax4, kde=False, bins=30, color='grey')
ax4.set_title('Pancreatic cancer dataset', fontsize=fs)
ax4.set_xlabel('Ratio of phylogenetically informative mutations to cells', fontsize=fs)
ax4.set_ylabel('Count', fontsize=fs)

# Plot for yang
sns.histplot(data=yang_data, 
            x='informative_muts_to_cell_ratio',
            ax=ax5, kde=False, bins=15, color='grey')
ax5.set_title('Lung cancer dataset (yang)', fontsize=fs)
ax5.set_xlabel('Ratio of phylogenetically informative mutations to cells', fontsize=fs)
ax5.set_ylabel('Count', fontsize=fs)

# Get the maximum y value across all plots
max_x = max(ax1.get_xlim()[1], ax2.get_xlim()[1], ax3.get_xlim()[1], ax4.get_xlim()[1], ax5.get_xlim()[1])

# Set the same y limits for all plots
ax1.set_xlim(0, max_x)
ax2.set_xlim(0, max_x)
ax3.set_xlim(0, max_x)
ax4.set_xlim(0, max_x)
ax5.set_xlim(0, max_x)

# ax3.set_ylim(0,10)
# ax3.set_yticks(range(0, 11, 1))

# Remove grid
ax1.grid(False)
ax2.grid(False)
ax3.grid(False)
ax4.grid(False)
ax5.grid(False)

ax1.tick_params(labelsize=fs)
ax2.tick_params(labelsize=fs)
ax3.tick_params(labelsize=fs)
ax4.tick_params(labelsize=fs)
ax5.tick_params(labelsize=fs)

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
