#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# Load the CSV file
file_path = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/model_selection_reseeding_no_reseeding_5_8_25_variable_rates_data_8_19_24/beam_classification_results.csv"

df = pd.read_csv(file_path)

# Classify beam Bayes factors as supporting reseeding or not
min_bf = 1.1  # 0 = support barely worth mentioning, 1.1 = positive support, 3 = strong support, 5 = overwhelming support
df["beam"] = df.apply(
    lambda row: (
        "yes" if row["bf"] >= min_bf
        else "no" if row["bf"] <= -min_bf
        else "nan"
    ),
    axis=1,
)

# Remove rows where beam is nan
df_original = df.copy()
df = df[df["beam"] != "nan"]

true_positive = np.sum((df["true"] == "yes") & (df["beam"] == "yes"))
true_negative = np.sum((df["true"] == "no") & (df["beam"] == "no"))
false_positive = np.sum((df["true"] == "no") & (df["beam"] == "yes"))
false_negative = np.sum((df["true"] == "yes") & (df["beam"] == "no"))

sensitivity = true_positive / (true_positive + false_negative)
specificity = true_negative / (true_negative + false_positive)

# do calculations with nan considered as no
false_negative = np.sum((df_original["true"] == "yes") & ((df_original["beam"] == "no") | (df_original["beam"] == "nan")))
sensitivity_nan = true_positive / (true_positive + false_negative)

# Plotting sensitivity and specificity bar plots
fig, ax = plt.subplots(figsize=(8, 6))
fs=18

ax.bar(['Sensitivity', 'Sensitivity\nwith intermediate\nas no-reseeding', 'Specificity'], [sensitivity, sensitivity_nan, specificity], color=['blue', 'orange', 'green'])
ax.set_ylim(0, 1)
ax.set_ylabel('Score', fontsize=fs)
ax.tick_params(axis="both", which="major", labelsize=fs)
plt.xticks(rotation=0, ha='center', fontsize=fs)

plt.tight_layout()
plt.savefig(file_path.replace(".csv", "_sensitivity_specificity.pdf"))
plt.close()
