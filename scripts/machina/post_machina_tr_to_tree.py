### Take in MACHINA internal node label output and leaf labeled tree and output tree with all node and leaf labels
from ete3 import Tree
import pandas as pd
import sys


def read_file_to_list_of_tuples(file_path):
    data = []
    node_num_to_name = {}
    node_num_to_color_code = {}
    color_code_to_tissue = {}
    with open(file_path, "r") as file:
        for line in file:
            # Split the line by spaces and create a tuple from the resulting elements
            elements = line.strip().split()
            if len(elements) == 2:
                name = elements[1].split('label="')[1]
                if "\\n" in name:
                    name = name.split("\\n")[0]
                else:
                    name = name.split('"')[0]
                node_num_to_name[elements[0]] = name

                color_code = elements[1].split("color=")[1].split(",")[0]
                node_num_to_color_code[elements[0]] = color_code

                if "\\n" in elements[1]:
                    tissue = elements[1].split("\\n")[1].split('"')[0]
                    color_code_to_tissue[color_code] = tissue
            elif len(elements) == 4:
                elements = (elements[0], elements[2])
                data.append(elements)
            else:
                continue

    edges = []
    for parent, child in data:
        parent_name = (
            node_num_to_name[parent]
            + "_"
            + color_code_to_tissue[node_num_to_color_code[parent]]
        )
        child_name = (
            node_num_to_name[child]
            + "_"
            + color_code_to_tissue[node_num_to_color_code[child]]
        )
        pair_renamed = (parent_name, child_name)
        edges.append(pair_renamed)
    return edges


leaf_tree = sys.argv[1]
machina_labels = sys.argv[2]
machina_dir = sys.argv[3]

# read in tree from edge list
connections = read_file_to_list_of_tuples(leaf_tree)
tree = Tree.from_parent_child_table(connections)

tree.write(format=8, outfile=f"{machina_dir}/machina_tree_all_tissue_labels.nwk")
