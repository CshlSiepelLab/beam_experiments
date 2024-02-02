### Take in MACHINA internal node label output and leaf labeled tree and output tree with all node and leaf labels
import ete3
import pandas as pd
import sys

leaf_tree = sys.argv[1]
machina_labels = sys.argv[2]
machina_dir = sys.argv[3]

# leaf_tree = "test_only_leaf_tissue_labels.nwk"
# machina_labels = "machina_results/T-t1-0.labeling"


tree = ete3.Tree(leaf_tree, format=8)
# Remove tissue labels for internal node names
for node in tree.traverse():
    if node.is_leaf() or node.is_root():
        continue
    else:
        current_name = node.name
        new_name = current_name.split("_")[0]
        node.name = new_name

tree.get_tree_root().name = '0'
machina_df = pd.read_csv(machina_labels, delim_whitespace = True, names = ['node', 'tissue'], index_col = 0)

for node in tree.traverse():
    if node.is_leaf() == False:
        node_name = int(node.name)
        tissue = machina_df.loc[node_name, 'tissue']
        new_name = str(node_name) + "_" + tissue
        node.name = new_name

tree.write(format=8, outfile = f'{machina_dir}/machina_tree_all_tissue_labels.nwk')