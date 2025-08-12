
import sys
import ete3
import pandas as pd
import random


def label_tissues_parsimony(tree, tissues_df, threshold_num_solutions):
    """
    Fitch parsimony algorithm to infer ancestral states of internal nodes for a tree
    """

    def postorder(node):
        if node.is_leaf():
            # Assign known tissue type from the tissues_df to the leaf node
            node.final_tissue = tissues_df.loc[
                tissues_df["cell"] == str(node.name), "tissue"
            ].values[0]
            node.tissue_set = {node.final_tissue}
            node.name = f"{node.name}_{node.final_tissue}"
            node.decision = "leaf"
        else:
            # Process all children
            children_tissue_sets = [postorder(child) for child in node.children]

            # Compute the possible tissues for internal nodes based on the children's tissue sets
            intersection = set.intersection(*children_tissue_sets)
            if intersection:
                node.tissue_set = intersection
            else:
                node.tissue_set = set.union(*children_tissue_sets)
        return node.tissue_set

    def preorder(node, total_solutions, parent_tissue=None):
        # Leaf tissues are already known so skip them for tissue assignment
        if not node.is_leaf():
            if node.is_root():
                # The root tissue is known
                node.final_tissue = f"{primary_tissue}"
                node.decision = "root"
            elif parent_tissue and parent_tissue in node.tissue_set:
                # If parent tissue is in the node's set, choose it
                node.final_tissue = parent_tissue
                node.decision = "parent"
            else:
                # If not then make an arbitrary choice from those available and increment the parsimony score
                num_tissues = len(list(node.tissue_set))
                if num_tissues == 1:
                    node.final_tissue = list(node.tissue_set)[0]
                    node.decision = "one_option"
                else:
                    node.final_tissue = random.choice(list(node.tissue_set))
                    node.decision = "random"
                total_solutions = total_solutions * num_tissues
                node.parsimony_score += 1
            node.name = f"{node.name}_{node.final_tissue}"
            # Recursively process children
            for child in node.children:
                total_solutions = preorder(child, total_solutions, node.final_tissue)
        else:
            # Check if leaf nodes are different tissues than their parents
            if parent_tissue != node.final_tissue:
                node.parsimony_score += 1
        return total_solutions

    def traverse_all_solutions(root):
        tree = root.copy()
        all_solutions = [tree]

        for node in tree.traverse():
            print(node.name)
            if node.decision == "random":
                print("random")
                tissues = list(node.tissue_set)
                split = node.name.split("_")
                label = split[0]
                first_tissue = split[1]
                tissues_subset = [tis for tis in tissues if tis != first_tissue]
                print(tissues_subset)
                count = len(all_solutions)
                for tissue in tissues_subset:
                    print(tissue)
                    for i in range(count):
                        print(i)
                        tree_copy = all_solutions[i].copy()
                        node_copy = tree_copy.search_nodes(name=node.name)[0]
                        node_copy.name = f"{label}_{tissue}"
                        all_solutions.append(tree_copy)

        return all_solutions

    # copy the input tree to avoid changing it in place
    copy_tree = tree.copy()

    # Run the postorder to get candidate tissues at each node
    postorder(copy_tree)

    # Initialize the parsimony scores
    for node in copy_tree.traverse():
        node.parsimony_score = 0

    # Assign the ancestral tissues for each node and update the parsimony score
    num_solutions = preorder(copy_tree, total_solutions=1)

    # Obtain the total parsimony score for the tree with random node selections
    total_parsimony_score = sum(node.parsimony_score for node in copy_tree.traverse())

    # If the total number of solutions is less than the specified threshold, then re-run the preorder and enumerate all solutions to be returned as a list of trees
    # print(f"Num solutions: {num_solutions}")
    if num_solutions < threshold_num_solutions and num_solutions != 1:
        all_solutions = traverse_all_solutions(copy_tree)
    elif num_solutions == 1:
        all_solutions = [copy_tree]
    else:
        all_solutions = []

    return copy_tree, all_solutions


# User inputs
tree_file = sys.argv[1]  # newick file
leaf_tissues_tsv = sys.argv[
    2
]  # tsv file with tip cell names and tissue labels as columns, and no header
outdir = sys.argv[3]  # where to write the output
primary_tissue = sys.argv[4]  # the known tissue label of the root node
threshold_num_solutions = int(
    sys.argv[5]
)  # the maximum number of possible solutions to enumerate them all in the output

tree = ete3.Tree(tree_file, format=8)
tissue_map = pd.read_csv(
    leaf_tissues_tsv,
    sep=r"\s+",
    header=None,
    names=["cell", "tissue"],
    dtype={"cell": str, "tissue": str},
)

# get results
random_parsimony_tree, all_solutions = label_tissues_parsimony(
    tree, tissue_map, threshold_num_solutions
)

parsimony_output = outdir + "/parsimony_tissues_random.nwk"
random_parsimony_tree.write(outfile=parsimony_output, format=8, format_root_node=True)

if len(all_solutions) > 0:
    # Write all solutions to a file
    for i, solution in enumerate(all_solutions):
        solution_output = f"{outdir}/parsimony_tissues_all_solutions_{i+1}.nwk"
        solution.write(outfile=solution_output, format=8, format_root_node=True)
