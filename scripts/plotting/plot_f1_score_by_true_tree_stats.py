
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


csv_path = sys.argv[1]
tree_stats_path = sys.argv[2]
output_dir = sys.argv[3]

# Read the CSV files into pandas DataFrames
df = pd.read_csv(csv_path)
tree_stats_df = pd.read_csv(tree_stats_path)

# Merge the dataframes on the 'sim_name' column
df = df.merge(tree_stats_df, left_on="sim", right_on="sim_name")

# Filter columns with "f1" in the name
f1_columns = [col for col in df.columns if "f1" in col and "PathFinder" not in col]

# Melt the dataframe to have 'f1' values in a single column
data_f1 = df.melt(
    id_vars=[
        "sim",
        "migration_count",
        "comigration_count",
        "num_multiedges",
        "met_to_met",
        "reseeding",
        "clonality",
    ],
    value_vars=f1_columns,
    var_name="method",
    value_name="f1",
)
data_f1["method"] = data_f1["method"].str.replace("_f1", "")
fs = 24

# Plot continuous variables
for column in ["migration_count", "comigration_count", "num_multiedges"]:
    g = sns.FacetGrid(data_f1, col="method", col_wrap=3, height=6, aspect=1)
    g.map_dataframe(sns.scatterplot, x=column, y="f1", alpha=0.5)
    g.map_dataframe(
        sns.kdeplot, x=column, y="f1", fill=True, alpha=0.3, common_norm=False
    )
    g.set_titles(col_template="{col_name}", size=fs)
    for ax in g.axes.flat:
        ax.tick_params(labelsize=fs)
    g.add_legend(title=None, fontsize=fs)
    g.set_axis_labels(column.replace("_", " ").title(), "F1 Score", fontsize=fs)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/f1_score_{column}.pdf")
    plt.close()

for column in ["migration_count", "comigration_count", "num_multiedges"]:
    data_f1[f"{column}_bin"] = pd.cut(data_f1[column], bins=5)
    data_f1[f"{column}_bin"] = data_f1[f"{column}_bin"].apply(
        lambda x: f"{int(x.left)}-{int(x.right)}" if pd.notnull(x) else "NaN"
    )
    plt.figure(figsize=(12, 8))
    sns.boxplot(x=f"{column}_bin", y="f1", hue="method", data=data_f1, showfliers=False)
    sns.stripplot(
        x=f"{column}_bin",
        y="f1",
        hue="method",
        data=data_f1,
        dodge=True,
        color="black",
        alpha=0.5,
        jitter=True,
        legend=False,
    )
    plt.ylabel("F1 Score", fontsize=fs)
    if column == "migration_count":
        x_label = "Migration count"
    elif column == "comigration_count":
        x_label = "Co-migration count"
    elif column == "num_multiedges":
        x_label = "Number of unique multi-edges"
    else:
        x_label = column.replace("_", " ").title()
    plt.xlabel(x_label, fontsize=fs)
    plt.xticks(rotation=0, fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(
        fontsize=fs,
        title=None,
        frameon=False,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    plt.tight_layout()
    plt.savefig(f"{output_dir}/f1_score_{column}_binned.pdf")
    plt.close()

# Plot binary variables
for column in ["met_to_met", "reseeding", "clonality"]:
    if column == "clonality":
        data_f1[column] = pd.Categorical(
            data_f1[column], categories=["Monoclonal", "Polyclonal"], ordered=True
        )
    else:
        data_f1[column] = pd.Categorical(
            data_f1[column], categories=[False, True], ordered=True
        )
    if column == "met_to_met":
        x_label = "Met to met"
    elif column == "reseeding":
        x_label = "Primary reseeding"
    elif column == "clonality":
        x_label = "Clonality"
    else:
        x_label = column.replace("_", " ").title()
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=column, y="f1", hue="method", data=data_f1, showfliers=False)
    sns.stripplot(
        x=column,
        y="f1",
        hue="method",
        data=data_f1,
        dodge=True,
        color="black",
        alpha=0.5,
        jitter=True,
        legend=False,
    )
    plt.ylabel("F1 Score", fontsize=fs)
    plt.xlabel(x_label, fontsize=fs)
    plt.xticks(rotation=0, fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(
        fontsize=fs,
        title=None,
        frameon=False,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    plt.tight_layout()
    plt.savefig(f"{output_dir}/f1_score_{column}.pdf")
    plt.close()
