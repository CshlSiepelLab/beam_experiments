
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt

serio_bf = "/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/marginal_likelihoods_gtr_random.csv"
serio_mi = "/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/gtr_beam_mutual_information.csv"
quinn_bf = "/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/marginal_likelihoods.csv"
quinn_mi = "/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/gtr_beam_mutual_information.csv"

outfile="/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/bf_mi_correlation.pdf"

# Read data
serio_bf_df = pd.read_csv(serio_bf)
serio_mi_df = pd.read_csv(serio_mi)
quinn_bf_df = pd.read_csv(quinn_bf)
quinn_mi_df = pd.read_csv(quinn_mi)

# Prepare dataframes with the same name for merging later
serio_bf_df["name"] = serio_bf_df["name"].str.replace("/", "_")
serio_mi_df["name"] = serio_mi_df["mouse_cp"].str.replace("/", "_")
quinn_bf_df["name"] = quinn_bf_df["name"].str.replace("/", "_")
quinn_mi_df["name"] = quinn_mi_df["mouse_cp"].str.replace("/", "_")

# Filter out columns
serio_mi_df = serio_mi_df[["name", "mutual_information_normalized"]]
serio_bf_df = serio_bf_df[["name", "bf(gtr-random)"]]
quinn_mi_df = quinn_mi_df[["name", "mutual_information_normalized"]]
quinn_bf_df = quinn_bf_df[["name", "bf(gtr-random)"]]

# Add 'dataset' column to each dataframe
serio_bf_df = serio_bf_df.rename(columns={"bf(gtr-random)": "serio_bf"})
serio_mi_df = serio_mi_df.rename(columns={"mutual_information_normalized": "serio_mi"})
quinn_bf_df = quinn_bf_df.rename(columns={"bf(gtr-random)": "quinn_bf"})
quinn_mi_df = quinn_mi_df.rename(columns={"mutual_information_normalized": "quinn_mi"})

# Merge dataframes
merged_df = serio_bf_df.merge(serio_mi_df, on="name", how="outer") \
    .merge(quinn_bf_df, on="name", how="outer") \
    .merge(quinn_mi_df, on="name", how="outer")
    
# Prepare data for each dataset
serio_df = merged_df.dropna(subset=["serio_bf", "serio_mi"])
quinn_df = merged_df.dropna(subset=["quinn_bf", "quinn_mi"])
all_df = pd.DataFrame({
    "Mutual Information Normalized": pd.concat([serio_df["serio_mi"], quinn_df["quinn_mi"]], ignore_index=True),
    "ln(Bayes factor)": pd.concat([serio_df["serio_bf"], quinn_df["quinn_bf"]], ignore_index=True),
    "Dataset": ["Serio"] * len(serio_df) + ["Quinn"] * len(quinn_df)
})

serio_df["Dataset"] = "Prostate cancer"
quinn_df["Dataset"] = "Lung cancer"

# Prepare merged dataframe for the right panel (no dataset distinction)
all_df = pd.concat([serio_df, quinn_df], ignore_index=True)
all_df = pd.DataFrame({
    "Mutual Information Normalized": pd.concat([serio_df["serio_mi"], quinn_df["quinn_mi"]], ignore_index=True),
    "ln(Bayes factor)": pd.concat([serio_df["serio_bf"], quinn_df["quinn_bf"]], ignore_index=True),
    "Dataset": ["Prostate cancer"] * len(serio_df) + ["Lung cancer"] * len(quinn_df)
})

fs = 20
fig, axes = plt.subplots(1, 3, figsize=(22, 6), sharex=False, sharey=False)

# Panel 1: Serio
sns.scatterplot(x=serio_df["serio_mi"], y=serio_df["serio_bf"], color="grey", ax=axes[0])
spearman_corr_serio, _ = spearmanr(serio_df["serio_mi"], serio_df["serio_bf"])
axes[0].set_xlabel("Mutual Information Normalized", fontsize=fs)
axes[0].set_ylabel("ln(Bayes factor)", fontsize=fs)
axes[0].set_title(f"Prostate cancer\nSpearman: {spearman_corr_serio:.2f}", fontsize=fs)
axes[0].tick_params(labelsize=fs)

# Panel 2: Quinn
sns.scatterplot(x=quinn_df["quinn_mi"], y=quinn_df["quinn_bf"], color="grey", ax=axes[1])
spearman_corr_quinn, _ = spearmanr(quinn_df["quinn_mi"], quinn_df["quinn_bf"])
axes[1].set_xlabel("Mutual Information Normalized", fontsize=fs)
axes[1].set_ylabel("ln(Bayes factor)", fontsize=fs)
axes[1].set_title(f"Lung cancer\nSpearman: {spearman_corr_quinn:.2f}", fontsize=fs)
axes[1].tick_params(labelsize=fs)

# Panel 3: Both datasets
sns.scatterplot(data=all_df, x="Mutual Information Normalized", y="ln(Bayes factor)", hue="Dataset", ax=axes[2], palette=["blue", "orange"])
spearman_corr_all, _ = spearmanr(all_df["Mutual Information Normalized"], all_df["ln(Bayes factor)"])
axes[2].set_xlabel("Mutual Information Normalized", fontsize=fs)
axes[2].set_ylabel("ln(Bayes factor)", fontsize=fs)
axes[2].set_title(f"Combined\nSpearman: {spearman_corr_all:.2f}", fontsize=fs)
axes[2].tick_params(labelsize=fs)
legend = axes[2].legend(fontsize=fs)
axes[2].legend(fontsize=fs, bbox_to_anchor=(1.15, 0.5), loc='center left', borderaxespad=0.)

plt.tight_layout()
plt.savefig(outfile)
plt.close()
