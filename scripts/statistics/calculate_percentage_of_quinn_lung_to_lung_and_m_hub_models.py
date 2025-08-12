
import sys
import csv
from matplotlib import pyplot as plt
import seaborn

consensus_graph_files = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/34/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/37/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/85/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/80/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/83/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/56/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/51/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/59/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/24/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/70/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/78/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/73/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/43/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/46/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/49/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/97/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/92/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/60/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/68/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/63/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/36/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/82/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/55/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/58/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/26/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/74/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/29/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/77/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/72/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/45/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/40/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/48/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/96/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/91/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/99/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/100/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/94/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/64/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/67/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/62/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/32/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/35/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/30/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/86/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/89/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/84/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/54/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/57/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/52/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/28/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/76/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/71/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/79/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/44/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/47/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/95/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/42/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/90/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/98/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/66/posterior_prob_graph.csv,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/61/posterior_prob_graph.csv"
outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/model_classifications_by_posterior_prob.csv"

primary_tissue = "LL"

files = consensus_graph_files.split(",")

thresholds = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

model_counts_by_threshold = {
    threshold: {"rl": 0, "m_only": 0, "others": 0, "none": 0}
    for threshold in thresholds
}

for file in files:

    if not open(file).read().strip():
        print(f"File is empty: {file}")
        continue

    for threshold in thresholds:

        tissues_seeded_from_primary = set()

        with open(file, "r") as f:
            for line in f.readlines():
                migration, prob = line.strip().split(",")
                prob = float(prob)
                source, target, num = migration.split("_")
                if prob >= threshold and source == primary_tissue:
                    tissues_seeded_from_primary.add(target)

        if "RL" in tissues_seeded_from_primary:
            model_counts_by_threshold[threshold]["rl"] += 1
        elif tissues_seeded_from_primary == {"M"}:
            model_counts_by_threshold[threshold]["m_only"] += 1
        elif tissues_seeded_from_primary:
            model_counts_by_threshold[threshold]["others"] += 1
        else:
            model_counts_by_threshold[threshold]["none"] += 1

# check the sums all match, so all solutions have been accounted for
total_files = len(files)
for threshold in thresholds:
    counts = list(model_counts_by_threshold[threshold].values())
    total_count = sum(counts)
    if total_count != total_files:
        print(f"Mismatch for threshold {threshold}: {total_count} != {total_files}")
        # write to output csv
        with open(outfile, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            csvfile.write("threshold,rl,m_only,others,none\n")
            for threshold in thresholds:
                writer.writerow(
                    [
                        threshold,
                        model_counts_by_threshold[threshold]["rl"],
                        model_counts_by_threshold[threshold]["m_only"],
                        model_counts_by_threshold[threshold]["others"],
                        model_counts_by_threshold[threshold]["none"],
                    ]
                )

# plot data for threshold groups as the x axis and the y axis as the counts with a single bar stacked for each group
thresholds = list(model_counts_by_threshold.keys())
rl = [model_counts_by_threshold[t]["rl"] for t in thresholds]
m_only = [model_counts_by_threshold[t]["m_only"] for t in thresholds]
others = [model_counts_by_threshold[t]["others"] for t in thresholds]
none = [model_counts_by_threshold[t]["none"] for t in thresholds]

# Create stacked bar plot
plt.figure(figsize=(12, 5))
fs = 22
x_positions = range(len(thresholds))  # Discrete x-axis positions

# Define colors for better aesthetics
colors = ["#4c72b0", "#dd8452", "#a3be8c", "darkgrey"]

# Create stacked bar plot with improved styling
plt.bar(x_positions, rl, color=colors[0], label="RL (and others)")
plt.bar(x_positions, m_only, bottom=rl, color=colors[1], label="M only")
plt.bar(
    x_positions,
    others,
    bottom=[i + j for i, j in zip(rl, m_only)],
    color=colors[2],
    label="Others",
)
plt.bar(
    x_positions,
    none,
    bottom=[i + j + k for i, j, k in zip(rl, m_only, others)],
    color=colors[3],
    label="None",
)

# Add labels, title, and legend with improved styling
plt.xlabel("Consensus Graph Threshold", fontsize=fs)
plt.ylabel("Percentage of All Data (%)", fontsize=fs)
plt.legend(
    fontsize=fs,
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False,
    title="Tissues seeded from primary",
    title_fontsize=fs,
)
plt.xticks(x_positions, thresholds, fontsize=fs)

# Convert y-axis to percentage
total_files = len(files)
y_ticks = range(0, 101, 25)  # Percentage ticks from 0 to 100
plt.yticks([y * total_files / 100 for y in y_ticks], y_ticks, fontsize=fs)

# Add gridlines for better readability
# plt.grid(axis="y", linestyle="--", alpha=0.7)

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig(outfile.replace(".csv", ".pdf"))
plt.close()
