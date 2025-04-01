#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


csv = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/pipelines/model_selection_mmus1467_cp01/beam_ns_mmus1467_cp01_no_reseeding_one_rate_reseeding_12_3_24/all_reps_results_ml_mmus1467_cp01.csv"

# Create DataFrame
df = pd.read_csv(csv)

# reformat model names
df["model"] = [model.replace("_", " ").capitalize() for model in df["model"]]

# Calculate means
means = df.groupby("model")["marginal_likelihood"].mean()
print("Means of each model:")
print(means)

# Calculate Bayes factor
bayes_factor = means["One rate reseeding"] - means["No reseeding"]
print(f"Bayes factor (One rate reseeding - No reseeding): {bayes_factor}")


# Plot boxplot
fs = 22

plt.figure(figsize=(6, 6))
sns.boxplot(
    x="model",
    y="marginal_likelihood",
    width=0.5,
    data=df,
    boxprops=dict(facecolor="none", edgecolor="black"),
    showfliers=False,
)
sns.stripplot(
    x="model",
    y="marginal_likelihood",
    data=df,
    jitter=True,
    size=5,
    color="grey",
    legend=False,
)
plt.ylabel("Marginal Likelihood", fontsize=fs)
plt.xlabel("")
plt.title(f"Mean Bayes factor = {bayes_factor:.2f}", fontsize=fs)
plt.xticks(fontsize=fs, rotation=45)
plt.yticks(fontsize=fs)

plt.tight_layout()
outfile = csv.replace(".csv", "_marginal_likelihoods_boxplot.pdf")
plt.savefig(outfile)
plt.close()
