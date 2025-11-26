import dendropy
import sys


def get_node_label(node):
    if node.taxon is not None:
        return node.taxon.label
    else:
        return node.label


def set_node_label(node, label):
    if node.taxon is not None:
        node.taxon.label = label
    else:
        node.label = label

tissue_tree_path = sys.argv[1]
branch_tree_path = sys.argv[2]
origin_tissue = sys.argv[3]
output_path = sys.argv[4]

tissue_tree = dendropy.Tree.get(path=tissue_tree_path, schema="newick", preserve_underscores=True)
branch_tree = dendropy.Tree.get(path=branch_tree_path, schema="newick", preserve_underscores=True)

branch_tree.suppress_unifurcations()    # Ensure no unifurcations due to laml formatting with root->origin

# Build tissue map from tissue tree
tissue_map = {}
for node in tissue_tree.preorder_node_iter():
    node_label = get_node_label(node)
    if "_" in node_label:
        base, tissue = node_label.rsplit("_", 1)
        tissue_map[base] = tissue

# Transfer tissue labels
for node in branch_tree.preorder_node_iter():
    if node.parent_node is None or node is branch_tree.seed_node:  # origin node
        set_node_label(node, f"_{origin_tissue}")
        continue
    node_label = get_node_label(node)
    if node_label in tissue_map:
        tissue = tissue_map[node_label]
        set_node_label(node, f"{node_label}_{tissue}")
    else:
        print("Warning: No tissue found for node", node_label)

# Write output newick
branch_tree.write(path=output_path, schema="newick", suppress_rooting=True, suppress_internal_node_labels=False)
