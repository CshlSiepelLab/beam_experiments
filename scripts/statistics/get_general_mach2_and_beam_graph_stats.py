import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy


mach2_file = "/grid/siepel/home/staklins/stored_results/beam/latest_results/general_graph_stats_for_beam_paper_from_latest_runs_8_2_25/mach2_all_results_8_18_25.csv"
beam_file = "/grid/siepel/home/staklins/stored_results/beam/latest_results/general_graph_stats_for_beam_paper_from_latest_runs_8_2_25/beam_all_results_8_18_25.csv"

# Output directory to save plots in
outdir = os.path.dirname(mach2_file)

# Read the input data in to dataframes
mach2_df = pd.read_csv(mach2_file)
beam_df = pd.read_csv(beam_file)

# Remove any CPs that were excluded in the paper, but may have slipped in here.
quinn_excluded_cps_5k_by_original_paper = ["5", "16", "18", "25", "33", "38", "39", "41", "50", "53", "65", "69", "75", "81", "87", "88", "93"]
quinn_excluded_cps_5k_by_compute_requirements = ["1", "2", "3", "4", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "19", "20", "21", "22", "23", "27", "31"]
all_quinn_excluded_cps = quinn_excluded_cps_5k_by_original_paper + quinn_excluded_cps_5k_by_compute_requirements
mach2_df = mach2_df[~((mach2_df["dataset_name"] == "quinn") & (mach2_df["cp"].isin(all_quinn_excluded_cps)))]
beam_df = beam_df[~((beam_df["dataset_name"] == "quinn") & (beam_df["cp"].isin(all_quinn_excluded_cps)))]

serio_cps_to_exlude = [("MMUS1875", "CP01"), ("MMUS1466", "CP01")]
mach2_df = mach2_df[~((mach2_df["dataset_name"] == "serio") & 
                        (mach2_df.apply(lambda row: (row["mouse"] , row["cp"]) in serio_cps_to_exlude, axis=1)))]
beam_df = beam_df[~((beam_df["dataset_name"] == "serio") & 
                    (beam_df.apply(lambda row: (row["mouse"] , row["cp"]) in serio_cps_to_exlude, axis=1)))]


# Remove CPs that were run with beam but not mach2, due to not having the primary tissue where mach2 requires it
mach2_mouse_cp = set(zip(mach2_df['mouse'], mach2_df['cp']))
beam_mouse_cp = set(zip(beam_df['mouse'], beam_df['cp']))
missing_in_mach2 = beam_mouse_cp - mach2_mouse_cp
beam_df = beam_df[~beam_df.apply(lambda row: (row['mouse'], row['cp']) in missing_in_mach2, axis=1)]

# Filter to only keep edges >= 0.5 for the comparison
mach2_df_05 = mach2_df[mach2_df["probability"] >= 0.5].copy()
beam_df_05 = beam_df[beam_df["probability"] >= 0.5].copy()

# Split 'source_target_edgenum' into 'source', 'target', and 'edgenum' cols
for df in [mach2_df_05, beam_df_05]:
    df[['source', 'target', 'edgenum']] = df['source_target_edgenum'].str.split('_', expand=True)
    df['edgenum'] = df['edgenum'].astype(int)

# Get migration counts
mach2_migration_counts = mach2_df_05.groupby(["dataset_name", "mouse", "cp"]).size().reset_index(name='migration_count')
beam_migration_counts = beam_df_05.groupby(["dataset_name", "mouse", "cp"]).size().reset_index(name='migration_count')
mach2_migration_counts_mean = mach2_migration_counts.groupby("dataset_name")["migration_count"].mean()
beam_migration_counts_mean = beam_migration_counts.groupby("dataset_name")["migration_count"].mean()

# Filter to only keep edges with edgenum == 1 to count co-migrations as unique edges in the graph
mach2_df_05_2 = mach2_df_05[mach2_df_05['edgenum'] == 1].copy()
beam_df_05_2 = beam_df_05[beam_df_05['edgenum'] == 1].copy()

# Get co-migration counts
mach2_comigration_counts = mach2_df_05_2.groupby(["dataset_name", "mouse", "cp"]).size().reset_index(name='comigration_count')
beam_comigration_counts = beam_df_05_2.groupby(["dataset_name", "mouse", "cp"]).size().reset_index(name='comigration_count')
mach2_comigration_counts_mean = mach2_comigration_counts.groupby("dataset_name")["comigration_count"].mean()
beam_comigration_counts_mean = beam_comigration_counts.groupby("dataset_name")["comigration_count"].mean()

# # Mean out degree (number of edges leaving each source tissue) and in degree (number of edges entering each target tissue)
# def mean_in_out_degree(df):
#     out_degree = df.groupby(["dataset_name", "mouse", "cp"])["source"].value_counts().groupby(["dataset_name", "mouse", "cp"]).mean().groupby("dataset_name").mean()
#     in_degree = df.groupby(["dataset_name", "mouse", "cp"])["target"].value_counts().groupby(["dataset_name", "mouse", "cp"]).mean().groupby("dataset_name").mean()
#     return out_degree, in_degree

# mach2_mean_out_degree, mach2_mean_in_degree = mean_in_out_degree(mach2_df_05)
# beam_mean_out_degree, beam_mean_in_degree = mean_in_out_degree(beam_df_05)

# Collect all stats into a DataFrame
rows = []
for method, migration_mean, comigration_mean in [
    ("mach2", mach2_migration_counts_mean, mach2_comigration_counts_mean),
    ("beam", beam_migration_counts_mean, beam_comigration_counts_mean)
]:
    for dataset in migration_mean.index:
        rows.append({
            "method": method,
            "dataset_name": dataset,
            "mean_migration_count": migration_mean[dataset],
            "mean_comigration_count": comigration_mean[dataset],
        })

stats_df = pd.DataFrame(rows)

# Format stats_df as a table and save as PDF
# Format column labels and cell text: replace _ with space and capitalize first letter
def format_label(label):
    return label.replace('_', ' ').capitalize()

def format_cell(cell):
    if isinstance(cell, str):
        return cell.replace('_', ' ').capitalize()
    return cell

# Some general formatting for the stats_df
stats_df["dataset_name"] = stats_df["dataset_name"].replace({"quinn": "lung cancer", "serio": "prostate cancer"})
stats_df = stats_df.sort_values(by=["dataset_name", "method"], key=lambda col: col if col.name != "method" else col.map({"mach2": 0, "beam": 1})).reset_index(drop=True)

# Format column labels
col_labels = [format_label(col) for col in stats_df.columns]

# Format cell text
cell_text = []
for row in stats_df.round(3).values:
    formatted_row = [format_cell(cell) for cell in row]
    cell_text.append(formatted_row)

fig, ax = plt.subplots(figsize=(10, 2 + 0.5 * len(stats_df)))
ax.axis('off')
tbl = ax.table(
    cellText=cell_text,
    colLabels=col_labels,
    loc='center',
    cellLoc='center'
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(10)
tbl.auto_set_column_width(col=list(range(len(stats_df.columns))))
plt.tight_layout()
pdf_path = os.path.join(outdir, "general_mach2_and_beam_graph_stats_table.pdf")
plt.savefig(pdf_path, bbox_inches='tight')
plt.close()
