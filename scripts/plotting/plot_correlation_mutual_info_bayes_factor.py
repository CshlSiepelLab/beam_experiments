
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt

# File paths
file1 = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/mi.csv"
file2 = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/bf.csv"
outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/correlation_mutual_info_bayes_factor.pdf"

# Read data
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Prepare and merge dataframes
df1["name"] = df1["mouse_cp"].str.replace("/", "_")
df2["name"] = df2["name"].str.replace("/", "_")

df1 = df1[["name", "mutual_information_normalized"]]
df2 = df2[["name", "bf(gtr-random)"]]

# Merge dataframes on 'sim_name'
merged_df = pd.merge(df1, df2, on="name")

# Extract relevant columns
mutual_info = merged_df["mutual_information_normalized"].astype(float)
bayes_factor = merged_df["bf(gtr-random)"].astype(float)

# Plotting
fs = 18

plt.figure(figsize=(6, 6))
sns.scatterplot(x=mutual_info, y=bayes_factor, color="grey")
# sns.regplot(x=mutual_info, y=bayes_factor, scatter=False, color='red')

# Calculate Spearman correlation
spearman_corr, _ = spearmanr(mutual_info, bayes_factor)

# # Linear regression
# X = mutual_info.values.reshape(-1, 1)
# y = bayes_factor.values
# reg = LinearRegression().fit(X, y)
# r_squared = reg.score(X, y)
# plt.plot(mutual_info, reg.predict(X), color='red', linewidth=2)


plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel("Mutual Information Normalized", fontsize=fs)
plt.ylabel("ln(Bayes factor)", fontsize=fs)
# plt.title(f'y = {reg.coef_[0]:.2f}x + {reg.intercept_:.2f}\nR^2: {r_squared:.2f}\nSpearman Correlation: {spearman_corr:.2f}', fontsize=fs)
plt.title(f"Spearman Correlation: {spearman_corr:.2f}", fontsize=fs)

plt.tight_layout()
plt.savefig(outfile)
plt.close()
