
import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import seaborn as sns

# user inputs
filepath=sys.argv[1]

# set outdir from input file
outdir = os.path.dirname(filepath)

# read in csv
df = pd.read_csv(filepath)
df["migration"] = df["dir_name"].str.split("/").apply(lambda x: x[-2])
df_long = pd.melt(
    df,
    id_vars=["dir_name", "migration"],
    value_vars=[
        "random_f1",
        "consensus_f1",
        "machina_f1",
        "beast_mcc_f1",
        "beast_posterior_f1",
    ],
)


# plot overall f1 score
textsize = 22
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
sns.boxplot(x="variable", y="value", data=df_long, width=0.5, color="lightblue")
plt.ylabel("F1 Score", fontsize=textsize)
plt.xlabel("")
plt.xticks(fontsize=textsize, rotation=-45)
plt.yticks(fontsize=textsize)
plt.grid(axis="y")
plt.tight_layout()
plt.savefig(f"{outdir}/f1_score_accuracy.pdf")

# plot migration f1 score
textsize = 22
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
sns.boxplot(
    x="variable",
    y="value",
    data=df_long,
    width=0.5,
    hue="migration",
    hue_order=["mS", "pS", "pM", "pR"],
)
plt.ylabel("F1 Score", fontsize=textsize)
plt.xlabel("")
plt.xticks(fontsize=textsize, rotation=-45)
plt.yticks(fontsize=textsize)
plt.grid(axis="y")
plt.tight_layout()
plt.savefig(f"{outdir}/migration_f1_score_accuracy.pdf")
