
import sys
import csv
from matplotlib import pyplot as plt
import seaborn

outfile = "/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/beam_gtr/5k/model_classifications_by_posterior_prob.csv"
infile_path = "/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/beam_gtr/5k"

cps=['34', '37', '80', '83', '56', '51', '59', '24', '70', '73', '43', '97', '92', '60', '68', '63', '36', '82', '55', '58', '26', '74', '77', '72', '45', '40', '91', '99', '100', '64', '67', '62', '32', '35', '30', '86', '89', '84', '54', '57', '52', '28', '76', '71', '79', '44', '47', '95', '42', '90', '98', '66', '61']

files=[]
for cp in cps:
    files.append(f"{infile_path}/{cp}/posterior_prob_graph.csv")

primary_tissue = "LL"

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

        if tissues_seeded_from_primary == {"RL", "M"} or tissues_seeded_from_primary == {"RL"}:
            model_counts_by_threshold[threshold]["rl"] += 1
        elif tissues_seeded_from_primary == {"M"}:
            model_counts_by_threshold[threshold]["m_only"] += 1
        elif tissues_seeded_from_primary:
            model_counts_by_threshold[threshold]["others"] += 1
        else:
            model_counts_by_threshold[threshold]["none"] += 1

# check the sums all match, so all solutions have been accounted for
with open(outfile, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["threshold", "rl", "m_only", "others", "none"])
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
plt.bar(x_positions, rl, color=colors[0], label="RL only or RL an M")
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
plt.xlabel("Posterior probability", fontsize=fs)
plt.ylabel("Percentage of CPs (%)", fontsize=fs)
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
