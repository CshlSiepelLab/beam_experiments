#!/usr/bin/env python3

import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import pickle


# Use to merge ideal and variable simulated rates data
ideal_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_2_20_25_uniform_50cells_50sites_data_7_24_24/precision_recall_curve/precision_recall_vars.pkl"
ideal_tree_stats_path = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_2_20_25_uniform_50cells_50sites_data_7_24_24/true_tree_stats.txt"

variable_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_2_25_25_data_from_8_19_24/precision_recall_curve/precision_recall_vars.pkl"
variable_tree_stats_path = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_2_25_25_data_from_8_19_24/true_tree_stats.txt"


outdir = os.path.dirname(ideal_file)

# Load data from pickle files
with open(ideal_file, "rb") as file:
    (
        ideal_machina_precisions,
        ideal_machina_recalls,
        ideal_metient_precisions,
        ideal_metient_recalls,
        ideal_random_precisions,
        ideal_random_recalls,
        ideal_consensus_precisions,
        ideal_consensus_recalls,
        ideal_parsimony_precisions,
        ideal_parsimony_recalls,
        ideal_fitchcount_precisions,
        ideal_fitchcount_recalls,
        ideal_mach2_all_thresh_df,
        ideal_all_thresh_df,
    ) = pickle.load(file)

with open(variable_file, "rb") as file:
    (
        variable_machina_precisions,
        variable_machina_recalls,
        variable_metient_precisions,
        variable_metient_recalls,
        variable_random_precisions,
        variable_random_recalls,
        variable_consensus_precisions,
        variable_consensus_recalls,
        variable_parsimony_precisions,
        variable_parsimony_recalls,
        variable_fitchcount_precisions,
        variable_fitchcount_recalls,
        variable_mach2_all_thresh_df,
        variable_all_thresh_df,
    ) = pickle.load(file)

# Concatenate data from ideal and variable inputs
machina_precisions = np.concatenate(
    (ideal_machina_precisions, variable_machina_precisions)
)
machina_recalls = np.concatenate((ideal_machina_recalls, variable_machina_recalls))
metient_precisions = np.concatenate(
    (ideal_metient_precisions, variable_metient_precisions)
)
metient_recalls = np.concatenate((ideal_metient_recalls, variable_metient_recalls))
random_precisions = np.concatenate(
    (ideal_random_precisions, variable_random_precisions)
)
random_recalls = np.concatenate((ideal_random_recalls, variable_random_recalls))
consensus_precisions = np.concatenate(
    (ideal_consensus_precisions, variable_consensus_precisions)
)
consensus_recalls = np.concatenate(
    (ideal_consensus_recalls, variable_consensus_recalls)
)
parsimony_precisions = np.concatenate(
    (ideal_parsimony_precisions, variable_parsimony_precisions)
)
parsimony_recalls = np.concatenate(
    (ideal_parsimony_recalls, variable_parsimony_recalls)
)
fitchcount_precisions = np.concatenate(
    (ideal_fitchcount_precisions, variable_fitchcount_precisions)
)
fitchcount_recalls = np.concatenate(
    (ideal_fitchcount_recalls, variable_fitchcount_recalls)
)
mach2_all_thresh_df = pd.concat(
    [ideal_mach2_all_thresh_df, variable_mach2_all_thresh_df]
)
all_thresh_df = pd.concat([ideal_all_thresh_df, variable_all_thresh_df])

# Load tree stats data
tree_stats_df = pd.concat(
    [pd.read_csv(ideal_tree_stats_path), pd.read_csv(variable_tree_stats_path)]
)

# Merge the dataframes on the 'sim_name' column
df = all_thresh_df.merge(tree_stats_df, left_on="sim", right_on="sim_name")
df2 = mach2_all_thresh_df.merge(tree_stats_df, left_on="sim", right_on="sim_name")

sim_names = all_thresh_df["sim"].unique()

fs = 24

# Plot continuous variables
for column in ["migration_count", "comigration_count", "num_multiedges"]:
    df[f"{column}_bin"] = pd.cut(df[column], bins=5)
    df[f"{column}_bin"] = df[f"{column}_bin"].apply(
        lambda x: f"{int(x.left)}-{int(x.right)}" if pd.notnull(x) else "NaN"
    )
    df2[f"{column}_bin"] = pd.cut(df2[column], bins=5)
    df2[f"{column}_bin"] = df2[f"{column}_bin"].apply(
        lambda x: f"{int(x.left)}-{int(x.right)}" if pd.notnull(x) else "NaN"
    )

    unique_bins = df[f"{column}_bin"].unique()
    num_bins = len(unique_bins)

    fig, axes = plt.subplots(1, num_bins, figsize=(5 * num_bins, 5), sharey=False)
    axes = axes.flatten()
    size = 100
    textsize = 22

    for i, bin_label in enumerate(
        sorted(unique_bins, key=lambda x: (int(x.split("-")[0]), int(x.split("-")[1])))
    ):
        ax = axes[i]
        bin_df = df[df[f"{column}_bin"] == bin_label]
        bin_df2 = df2[df2[f"{column}_bin"] == bin_label]
        bin_sims = bin_df["sim"].unique()
        num_sims = len(bin_sims)

        avg_df = (
            bin_df.groupby("Threshold")[["precision", "recall"]].mean().reset_index()
        )
        avg_df2 = (
            bin_df2.groupby("Threshold")[["precision", "recall"]].mean().reset_index()
        )

        avg_machina_precision = np.nanmean(
            [machina_precisions[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_machina_recall = np.nanmean(
            [machina_recalls[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_metient_precision = np.nanmean(
            [metient_precisions[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_metient_recall = np.nanmean(
            [metient_recalls[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_random_precision = np.nanmean(
            [random_precisions[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_random_recall = np.nanmean(
            [random_recalls[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_consensus_precision = np.nanmean(
            [consensus_precisions[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_consensus_recall = np.nanmean(
            [consensus_recalls[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_parsimony_precision = np.nanmean(
            [parsimony_precisions[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_parsimony_recall = np.nanmean(
            [parsimony_recalls[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_fitchcount_precision = np.nanmean(
            [fitchcount_precisions[sim_names.tolist().index(sim)] for sim in bin_sims]
        )
        avg_fitchcount_recall = np.nanmean(
            [fitchcount_recalls[sim_names.tolist().index(sim)] for sim in bin_sims]
        )

        if not avg_df.empty:
            ax.plot(avg_df["recall"], avg_df["precision"], color="grey", label="BEAM")
        if not avg_df2.empty:
            ax.plot(
                avg_df2["recall"], avg_df2["precision"], color="navy", label="MACH2"
            )
        if not np.isnan(avg_machina_recall) and not np.isnan(avg_machina_precision):
            ax.scatter(
                avg_machina_recall,
                avg_machina_precision,
                color="red",
                label="MACHINA",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_metient_recall) and not np.isnan(avg_metient_precision):
            ax.scatter(
                avg_metient_recall,
                avg_metient_precision,
                color="green",
                label="Metient",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_consensus_recall) and not np.isnan(avg_consensus_precision):
            ax.scatter(
                avg_consensus_recall,
                avg_consensus_precision,
                color="blue",
                label="Consensus",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_random_recall) and not np.isnan(avg_random_precision):
            ax.scatter(
                avg_random_recall,
                avg_random_precision,
                color="black",
                label="Random",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_parsimony_recall) and not np.isnan(avg_parsimony_precision):
            ax.scatter(
                avg_parsimony_recall,
                avg_parsimony_precision,
                color="purple",
                label="Parsimony",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_fitchcount_recall) and not np.isnan(
            avg_fitchcount_precision
        ):
            ax.scatter(
                avg_fitchcount_recall,
                avg_fitchcount_precision,
                color="orange",
                label="FitchCount",
                s=size,
                marker="x",
            )

        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel("Recall", fontsize=textsize)
        ax.set_ylabel("Precision", fontsize=textsize)
        ax.set_title(f"{bin_label}\n(n={num_sims})", fontsize=textsize)
        ax.tick_params(axis="both", which="major", labelsize=textsize)
        if i == num_bins - 1:
            ax.legend(
                bbox_to_anchor=(1.05, 0.5),
                loc="upper left",
                fontsize=14,
                edgecolor="none",
            )
    if column == "migration_count":
        title = "Migration count"
    elif column == "comigration_count":
        title = "Co-migration count"
    elif column == "num_multiedges":
        title = "Number of unique multi-edges"
    fig.suptitle(title, fontsize=fs)
    plt.tight_layout()
    outfile = (
        f"{outdir}/precision_recall_by_{column}_bin_merged_all_ideal_variable_sims.pdf"
    )
    plt.savefig(outfile, bbox_inches="tight")
    plt.close()


# Plot binary variables
for column in ["met_to_met", "reseeding", "clonality"]:
    if column == "clonality":
        unique_values = ["Monoclonal", "Polyclonal"]
    else:
        unique_values = [False, True]
    num_values = len(unique_values)

    fig, axes = plt.subplots(1, num_values, figsize=(5 * num_values, 5), sharey=False)
    axes = axes.flatten()
    size = 100
    textsize = 22

    for i, value in enumerate(unique_values):
        ax = axes[i]
        value_df = df[df[column] == value]
        value_df2 = df2[df2[column] == value]
        value_sims = value_df["sim"].unique()
        num_sims = len(value_sims)

        avg_df = (
            value_df.groupby("Threshold")[["precision", "recall"]].mean().reset_index()
        )
        avg_df2 = (
            value_df2.groupby("Threshold")[["precision", "recall"]].mean().reset_index()
        )

        avg_machina_precision = np.nanmean(
            [machina_precisions[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_machina_recall = np.nanmean(
            [machina_recalls[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_metient_precision = np.nanmean(
            [metient_precisions[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_metient_recall = np.nanmean(
            [metient_recalls[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_random_precision = np.nanmean(
            [random_precisions[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_random_recall = np.nanmean(
            [random_recalls[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_consensus_precision = np.nanmean(
            [consensus_precisions[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_consensus_recall = np.nanmean(
            [consensus_recalls[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_parsimony_precision = np.nanmean(
            [parsimony_precisions[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_parsimony_recall = np.nanmean(
            [parsimony_recalls[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_fitchcount_precision = np.nanmean(
            [fitchcount_precisions[sim_names.tolist().index(sim)] for sim in value_sims]
        )
        avg_fitchcount_recall = np.nanmean(
            [fitchcount_recalls[sim_names.tolist().index(sim)] for sim in value_sims]
        )

        if not avg_df.empty:
            ax.plot(avg_df["recall"], avg_df["precision"], color="grey", label="BEAM")
        if not avg_df2.empty:
            ax.plot(
                avg_df2["recall"], avg_df2["precision"], color="navy", label="MACH2"
            )
        if not np.isnan(avg_machina_recall) and not np.isnan(avg_machina_precision):
            ax.scatter(
                avg_machina_recall,
                avg_machina_precision,
                color="red",
                label="MACHINA",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_metient_recall) and not np.isnan(avg_metient_precision):
            ax.scatter(
                avg_metient_recall,
                avg_metient_precision,
                color="green",
                label="Metient",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_consensus_recall) and not np.isnan(avg_consensus_precision):
            ax.scatter(
                avg_consensus_recall,
                avg_consensus_precision,
                color="blue",
                label="Consensus",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_random_recall) and not np.isnan(avg_random_precision):
            ax.scatter(
                avg_random_recall,
                avg_random_precision,
                color="black",
                label="Random",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_parsimony_recall) and not np.isnan(avg_parsimony_precision):
            ax.scatter(
                avg_parsimony_recall,
                avg_parsimony_precision,
                color="purple",
                label="Parsimony",
                s=size,
                marker="x",
            )
        if not np.isnan(avg_fitchcount_recall) and not np.isnan(
            avg_fitchcount_precision
        ):
            ax.scatter(
                avg_fitchcount_recall,
                avg_fitchcount_precision,
                color="orange",
                label="FitchCount",
                s=size,
                marker="x",
            )

        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel("Recall", fontsize=textsize)
        ax.set_ylabel("Precision", fontsize=textsize)
        ax.set_title(f"{value}\n(n={num_sims})", fontsize=textsize)
        ax.tick_params(axis="both", which="major", labelsize=textsize)
        if i == num_values - 1:
            ax.legend(
                bbox_to_anchor=(1.05, 0.5),
                loc="upper left",
                fontsize=14,
                edgecolor="none",
            )

    if column == "met_to_met":
        title = "Met to met"
    elif column == "reseeding":
        title = "Primary reseeding"
    elif column == "clonality":
        title = "Clonality"
    fig.suptitle(title, fontsize=fs)
    plt.tight_layout()
    outfile = (
        f"{outdir}/precision_recall_by_{column}_merged_all_ideal_variable_sims.pdf"
    )
    plt.savefig(outfile, bbox_inches="tight")
    plt.close()
