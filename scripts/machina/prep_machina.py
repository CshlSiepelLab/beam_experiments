### Given a newick file with leaves labeled by name_tissue, then prep input Tree edge and leaf label files for MACHINA
import sys
import ete3
import pandas as pd

leaf_labeled_tree = sys.argv[1]
output_dir = sys.argv[2]

# leaf_labeled_tree = "test_only_leaf_tissue_labels.nwk"

input_file = leaf_labeled_tree.split("/")[-1]
input_prefix = input_file.split(".")[0]
output_file_leaf = output_dir + "/" + input_prefix + ".labeling"
output_file_edges = output_dir + "/" + input_prefix + ".tree"
output_file_colors = output_dir + "/" + input_prefix + "_colors.txt"

tree = ete3.Tree(leaf_labeled_tree, format=8)
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

tissues = leaf_label['tissue'].unique()
i = 1
color_map = {}
for tissue in tissues:
    color_map[tissue] = i
    i += 1

color_map_data = {'tissue': tissues, 'color': range(len(tissues))}
color_map = pd.DataFrame(color_map_data)

leaf_label.to_csv(output_file_leaf, sep="\t", index=False, header = False)
edges.to_csv(output_file_edges, sep="\t", index=False, header = False)
color_map.to_csv(output_file_colors, sep="\t", index=False, header = False)

