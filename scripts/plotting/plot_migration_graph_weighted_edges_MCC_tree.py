
import re, sys
import random
import pandas as pd
import numpy as np
import networkx as nx
import ete3
import matplotlib.pyplot as plt
import matplotlib


def remove_bracket_content(match):
    annotations.append(match.group()[1:-1])
    return ""


def get_labels(newick):
    label_pattern = re.compile(r"([A-Za-z0-9_]+):")
    leaf_labels = label_pattern.findall(newick)
    return leaf_labels


def label_nodes(newick):
    leaf_labels = get_labels(newick)
    try:
        leaf_labels_max = max(map(int, leaf_labels))
    except ValueError:
        # Handle the case where leaf labels are not convertible to integers
        leaf_labels_max = len(leaf_labels)
    # add labels to nodes
    start_label = leaf_labels_max + 1
    node_labeled_newick = ""
    parts = newick.split(")")
    for part in parts[:-1]:
        part = part + ")" + str(start_label)
        node_labeled_newick += part
        start_label += 1
    node_labeled_newick += parts[-1]
    return node_labeled_newick


consensus_tree_file = sys.argv[1]
primary_tissue = sys.argv[2]


with open(consensus_tree_file, "r") as file:
    for line in file:
        line = line.strip()
        if line.startswith("tree"):
            tree_info = line

# remove tree name and = from tree info to get only the newick
tree = "".join(tree_info.split(" ")[3:])

# strip tree string to newick with associated dataframe of annotations
bracket_content_pattern = re.compile(r"\[.*?\]")
annotations = []
newick = re.sub(bracket_content_pattern, remove_bracket_content, tree)
annotations = [re.split(r",(?![^{]*})", x.replace("&", "")) for x in annotations]

annotations = [
    {
        key: value.replace("{", "").replace("}", "")
        for trait in annotation
        for key, value in [trait.split("=")]
    }
    for annotation in annotations
]

# label nodes in newick with only leaf labels
node_labeled_newick = label_nodes(newick)
node_labels = get_labels(node_labeled_newick)

# make a dictionary for annotations to node labels
annotations_dict = {}
for node in node_labels:
    annotations_dict[node] = annotations[node_labels.index(node)]

annotations_df = pd.DataFrame.from_dict(annotations_dict, orient="index")
location_probs_df = annotations_df.loc[:, ["location.set", "location.set.prob"]]

# read newick into ete3 Tree
tree = ete3.Tree(node_labeled_newick, format=3)

locations = list(annotations_df["location.set"].values)
uniq_locations = list(
    set([value for location in locations for value in location.split(",")])
)

# make adjacency matrix weighted by location probabilities with source as row index names and recipient as column names
weighted_adjacency_matrix = pd.DataFrame(
    0, index=uniq_locations, columns=uniq_locations
)

for node in tree.traverse():
    if node.is_root():
        continue
    node_name = node.name
    parent_name = node.up.name
    node_locs = list(location_probs_df.loc[node_name, "location.set"].split(","))
    node_locs_probs = list(
        location_probs_df.loc[node_name, "location.set.prob"].split(",")
    )
    parent_locs = list(location_probs_df.loc[parent_name, "location.set"].split(","))
    parent_locs_probs = list(
        location_probs_df.loc[parent_name, "location.set.prob"].split(",")
    )
    for parent_loc in parent_locs:
        parent_loc_index = parent_locs.index(parent_loc)
        parent_loc_prob = float(parent_locs_probs[parent_loc_index])
        for node_loc in node_locs:
            node_loc_index = node_locs.index(node_loc)
            node_loc_prob = float(node_locs_probs[node_loc_index])
            joint_prob = parent_loc_prob * node_loc_prob
            weighted_adjacency_matrix.loc[parent_loc, node_loc] = (
                weighted_adjacency_matrix.loc[parent_loc, node_loc] + joint_prob
            )

# remove diagonal entries to not plot self-migrations
num_tissues = len(weighted_adjacency_matrix)
for i in range(0, num_tissues):
    for j in range(0, num_tissues):
        if i == j:
            weighted_adjacency_matrix.iloc[i, j] = 0

# normalize weighted adjacency matrix so the largest is 1
max_value = weighted_adjacency_matrix.values.max()
weighted_adj_norm = weighted_adjacency_matrix / max_value

# make complete graph with all locations
H = nx.MultiDiGraph()
for loc1 in uniq_locations:
    for loc2 in uniq_locations:
        if loc1 != loc2:
            weight = weighted_adj_norm.loc[loc1, loc2]
            H.add_edge(loc1, loc2, weight=weight)

nodes_np = sorted(list(H.nodes()))
nodes = [primary_tissue] + [node for node in nodes_np if node != primary_tissue]

G = nx.MultiDiGraph()
G.add_nodes_from(nodes)
G.add_edges_from(H.edges(data=True))

# Extract edge weights
edge_weights = [G.get_edge_data(u, v)[0]["weight"] for u, v in G.edges()]
edge_widths = [weight * 5 for weight in edge_weights]

# Create a colormap based on edge weights
# cmap = matplotlib.colormaps['binary']
cmap = plt.cm.colors.LinearSegmentedColormap.from_list(
    "custom", [(1, 1, 1), (1, 1, 1), (0.5, 0.5, 0.5), (0, 0, 0)], N=256
)

# Draw the graph with edge colors
edge_colors = [sm.to_rgba(weight) for (u, v), weight in zip(G.edges(), edge_weights)]
node_colors = range(len(nodes))
node_cmap = matplotlib.cm.get_cmap("tab20", len(nodes))

# Find the node corresponding to the primary tissue
primary_tissue_node = [node for node in G.nodes() if node == primary_tissue][0]

# Create positions for the nodes
fig, ax = plt.subplots(figsize=(8, 8))
max_width = ax.get_position().width
pos = {}
row_height = 0.1
num_nodes = len(nodes)

for i, node in enumerate(G.nodes()):
    if node == primary_tissue:
        pos[node] = (max_width / 2, 0)
    else:
        pos[node] = (
            (max_width / num_nodes) * (i + 0.5),
            -row_height + random.uniform(0, 0.025),
        )

# make my own color map of 10 colors for now
node_colors = [
    "black",
    "red",
    "green",
    "blue",
    "orange",
    "purple",
    "brown",
    "pink",
    "gray",
    "gold",
]
node_colors = node_colors[0 : len(nodes)]

nx.draw(
    G,
    pos=pos,
    ax=ax,
    with_labels=False,
    connectionstyle="arc3, rad = 0.2",
    edge_color=edge_colors,
    edge_cmap=cmap,
    width=edge_widths,
    arrowsize=20,
    font_size=10,
    font_color="black",
    font_weight="bold",
    node_shape="s",
    node_size=1000,
    node_color=node_colors,
)
# cmap=node_cmap)

# Add a colorbar to show the weight gradient
sm = plt.cm.ScalarMappable(cmap=cmap)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.5, ticks=np.arange(0, 1.01, 0.25))
cbar.set_label("Probability")

# Create legend for node colors
# legend_labels = {loc: node_cmap(i) for i, loc in enumerate(list(G.nodes()))}
legend_labels = {loc: node_color for loc, node_color in zip(nodes, node_colors)}
legend_handles = [
    plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=color, markersize=10)
    for color in legend_labels.values()
]
ax.legend(
    legend_handles,
    legend_labels.keys(),
    title="Node Locations",
    loc="upper left",
    bbox_to_anchor=(0.9, 1),
)


# plt.show()

outfile = consensus_tree_file.split(".")[0] + "_mcc_migration_graph.pdf"
plt.savefig(outfile)
