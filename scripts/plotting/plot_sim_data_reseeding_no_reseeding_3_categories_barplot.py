import matplotlib.pyplot as plt
import numpy as np

# Data preparation
categories = ['No PR\n(lBf < 3)', 'Between\nclassification\nthresholds', 'PR\n(lBf > 3)']

# Data for each category
no_counts = [5, 34, 0]
yes_counts = [0, 29, 10]

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Width of bars
bar_width = 0.5
x_pos = np.arange(len(categories))

# Create stacked bar chart
bars_no = ax.bar(x_pos, no_counts, bar_width, label='No PR', color='orange', alpha=0.8)
bars_yes = ax.bar(x_pos, yes_counts, bar_width, bottom=no_counts, label='PR', color='blue', alpha=0.8)

# Add value labels on bars
grand_total = sum(no_counts) + sum(yes_counts)

for i, (no_count, yes_count) in enumerate(zip(no_counts, yes_counts)):
    total = no_count + yes_count

    # Calculate percentages
    no_pct_within = 100 * no_count / total if total > 0 else 0
    yes_pct_within = 100 * yes_count / total if total > 0 else 0

    # Label for 'No' portion (if present)
    if no_count > 0:
        ax.text(i, no_count/2, f'{no_count}\n({no_pct_within:.1f}%)', 
                ha='center', va='center', fontsize=14)

    # Label for 'Yes' portion (if present)
    if yes_count > 0:
        ax.text(i, no_count + yes_count/2, f'{yes_count}\n({yes_pct_within:.1f}%)', 
                ha='center', va='center', fontsize=14)

    # Total label above bar
    ax.text(i, total + 1, f'Total: {total}\n({100*total/grand_total:.1f}%)', 
            ha='center', va='bottom', fontsize=14)

ax.legend(title="True PR classification", fontsize=18, title_fontsize=18, frameon=False)
ax.set_xlabel('', fontsize=18)
ax.set_ylabel('Count', fontsize=18)
ax.set_xticks(x_pos)
ax.set_xticklabels(categories, fontsize=18)
ax.tick_params(axis='y', labelsize=18)

# Set y-axis limit to accommodate labels
ax.set_ylim(0, max(np.array(no_counts) + np.array(yes_counts)) + 10)

# Adjust layout to prevent clipping
plt.tight_layout()

# Display the plot
plt.savefig('/grid/siepel/home/staklins/stored_results/bayesian_migration_graph_inference/latest_results/model_selection_reseeding_no_reseeding_5_8_25_variable_rates_data_8_19_24/3_category_barplot.pdf', dpi=300, bbox_inches='tight')
plt.close()
