import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

df = pd.read_csv("/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_24_24/precision_recall_curve/metient_all_threshold_stats.csv")



# Filter thresholds between 0 and 1 in 0.1 increments
thresholds = np.arange(0, 1.01, 0.1).round(2).tolist()
df['Threshold'] = df['Threshold'].astype(float).round(2)
# Round both sides for safe matching
df_subset = df[df['Threshold'].isin(thresholds)]

# Initialize figure
plt.figure(figsize=(12, 5))

# Alternatively, if you want both precision and recall side-by-side per threshold:
df_melted = df_subset.melt(id_vars="Threshold", value_vars=["precision", "recall"],
                           var_name="Metric", value_name="Score")
sns.boxplot(data=df_melted, x="Threshold", y="Score", hue="Metric", fliersize=0)
sns.stripplot(data=df_melted, x="Threshold", y="Score", hue="Metric",
              dodge=True, alpha=0.4, jitter=0.2, marker='o', linewidth=0.5, size=3, legend=False)

# Add mean points
group_means = df_melted.groupby(["Threshold", "Metric"], as_index=False)["Score"].mean()
sns.pointplot(data=group_means, x="Threshold", y="Score", hue="Metric",
              dodge=0.4, markers="D", scale=0.8, errwidth=0, linestyles="", color="black")


plt.xlabel("Threshold")
plt.ylabel("Score")
plt.title("Precision and Recall Distributions Across Thresholds")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend(bbox_to_anchor=(1.02, 0.5))
plt.ylim(0,1)
plt.tight_layout()

# Save instead of showing
plt.savefig("/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_24_24/precision_recall_curve/explore_metient_curve_vs_point_discrepancy_11_3_25/metient_at_diff_thresholds.pdf")
plt.close()

# Plot best solution precision/recall
best_df = pd.read_csv("/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_24_24/precision_recall_curve/metrics.csv")
best_df = best_df[["Metient_bestSol_precision", "Metient_bestSol_recall"]]

best_df = best_df.melt(var_name="Metric", value_name="Score")

plt.figure(figsize=(6, 6))
sns.boxplot(data=best_df, x="Metric", y="Score", fliersize=0, palette="pastel")
sns.stripplot(data=best_df, x="Metric", y="Score", alpha=0.4, jitter=0.2, marker='o', linewidth=0.5, size=4)

# Add mean points
means = best_df.groupby("Metric", as_index=False)["Score"].mean()
sns.pointplot(data=means, x="Metric", y="Score", markers="D", color="black", scale=1.0, errwidth=0, linestyles="")


plt.ylabel("Score")
plt.title("Metient Best-Solution Precision and Recall")
plt.ylim(0, 1)
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_24_24/precision_recall_curve/explore_metient_curve_vs_point_discrepancy_11_3_25/metient_best_solutions.pdf")
plt.close()
