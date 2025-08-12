
import sys
import ete3
import pandas as pd

leaf_labeled_tree = sys.argv[1]
output_dir = sys.argv[2]
primary_tissue = str(sys.argv[3])
leaf_labels_tsv = str(sys.argv[4])

input_file = leaf_labeled_tree.split("/")[-1]
input_prefix = input_file.split(".")[0]
output_file_leaf = output_dir + "/" + input_prefix + ".labeling"
output_file_edges = output_dir + "/" + input_prefix + ".tree"
output_file_colors = output_dir + "/" + input_prefix + "_colors.txt"

tree = ete3.Tree(leaf_labeled_tree, format=8)

# Remove tissue labels for internal node names, if they exist and add names to nodes that are dummy nodes for identical sequences from laml output
i=1
for node in tree.traverse():
    if node.is_root():
        node.name = "root"
    elif not node.is_leaf():
        current_name = node.name
        new_name = current_name.split("_")[0]
        if new_name == "node" and current_name.split("_")[1] == "dummy":
            new_name = f"nodedummy{i}"
            i += 1
        node.name = new_name

edges = []

for node in tree.traverse():
    if not node.is_leaf():
        for child in node.children:
            edges.append((node.name, child.name))

leaf_label = pd.read_csv(leaf_labels_tsv, sep="\s+", names=["leaf", "tissue"])
tissues = leaf_label["tissue"].unique().tolist()

# Fix when primary tissue is not a leaf label, but required in coloring scheme for MACHINA to run
if primary_tissue not in tissues:
    tissues.append(primary_tissue)

i = 1
color_map = {}
for tissue in tissues:
    color_map[tissue] = i
    i += 1

# output files
with open(output_file_edges, "w") as file:
    for edge in edges:
        file.write(f"{edge[0]}\t{edge[1]}\n")

leaf_label.to_csv(output_file_leaf, sep="\t", index=False, header=False)

with open(output_file_colors, "w") as file:
    for key, value in color_map.items():
        file.write(f"{key}\t{value}\n")
