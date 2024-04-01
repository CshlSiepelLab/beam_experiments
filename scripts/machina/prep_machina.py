### Given a newick file with leaves labeled by name_tissue, then prep input Tree edge and leaf label files for MACHINA
import sys
import ete3
import pandas as pd

leaf_labeled_tree = sys.argv[1]
output_dir = sys.argv[2]
primary_tissue = str(sys.argv[3])

# leaf_labeled_tree = "machina_m8_sim_data/seed10046/T_seed10046_tissue_labeled_true_tree.nwk"
# output_dir = "machina_m8_sim_data/seed10046/machina"
# primary_tissue = "P"

input_file = leaf_labeled_tree.split("/")[-1]
input_prefix = input_file.split(".")[0]
output_file_leaf = output_dir + "/" + input_prefix + ".labeling"
output_file_edges = output_dir + "/" + input_prefix + ".tree"
output_file_colors = output_dir + "/" + input_prefix + "_colors.txt"

tree = ete3.Tree(leaf_labeled_tree, format=8)


# Remove tissue labels for internal node names
for node in tree.traverse():
    if node.is_leaf() or node.is_root():
        continue
    else:
        current_name = node.name
        new_name = current_name.split("_")[0]
        node.name = new_name

tree.get_tree_root().name = '0'


leaf_label = pd.DataFrame(columns = ['leaf', 'tissue'])
edges = pd.DataFrame(columns = ['node1', 'node2'])

for node in tree.traverse():
    if node.is_leaf() == True:
        name,tissue = node.name.split("_")
        leaf_label.loc[len(leaf_label)] = [name, tissue]
    else:
        node_name = node.name
        children = node.children
        for child in children:
            child_name = child.name
            if "_" in child_name:
                child_name = child_name.split("_")[0]
            edges.loc[len(edges)] = [node_name, child_name]

tissues = leaf_label['tissue'].unique().tolist()
# Fix when primary tissue is not a leaf label, but required in coloring scheme for MACHINA to run
if primary_tissue not in tissues:
    tissues.append(primary_tissue)
i = 1
color_map = {}
for tissue in tissues:
    color_map[tissue] = i
    i += 1

# Output files for MACHINA
leaf_label.to_csv(output_file_leaf, sep="\t", index=False, header = False)
edges.to_csv(output_file_edges, sep="\t", index=False, header = False)

with open(output_file_colors, "w") as file:
    for key,value in color_map.items():
        file.write(f'{key}\t{value}\n')


