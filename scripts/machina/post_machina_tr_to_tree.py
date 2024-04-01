### Take in MACHINA internal node label output and leaf labeled tree and output tree with all node and leaf labels
from ete3 import Tree
import pandas as pd
import sys

def read_file_to_list_of_tuples(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line by spaces and create a tuple from the resulting elements
            elements = tuple(line.strip().split())
            data.append(elements)
    return data

leaf_tree = sys.argv[1]
machina_labels = sys.argv[2]
machina_dir = sys.argv[3]

# leaf_tree = "machina_m5_sim_data/seed955/machina/T-P-0.tree"
# machina_labels = "machina_m5_sim_data/seed955/machina/T-P-0.labeling"
# machina_dir = "machina_m5_sim_data/seed955/machina"

# read in tree from edge list
connections = read_file_to_list_of_tuples(leaf_tree)
tree = Tree.from_parent_child_table(connections)

# Remove tissue labels for internal node names
for node in tree.traverse():
    if node.is_leaf() or node.is_root():
        continue
    else:
        current_name = node.name
        new_name = current_name.split("_")[0]
        node.name = new_name

tree.get_tree_root().name = '0'
machina_df = pd.read_csv(machina_labels, delim_whitespace = True, names = ['node', 'tissue'], dtype={'node':str,'tissue':str})

for node in tree.traverse():
    node_name = node.name
    row = machina_df.loc[machina_df['node'] == node_name]
    tissue = row['tissue'].values[0]
    new_name = str(node_name) + "_" + tissue
    node.name = new_name

tree.write(format=8, outfile = f'{machina_dir}/machina_tree_all_tissue_labels.nwk')