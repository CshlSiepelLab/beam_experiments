#!/usr/bin/env python3

import sys
import os
from ete3 import Tree, TreeStyle, NodeStyle, CircleFace, TextFace

def extract_tissue_from_name(tip_name):
    """
    Extract tissue from tip name in format 'tissue.name'
    Returns tissue and name as tuple
    """
    if '.' in tip_name:
        parts = tip_name.split('.', 1)  # Split on first '.' only
        if len(parts) == 2:
            tissue = parts[0]
            # Remove numbers from tissue name for coarse grained annotations
            tissue = ''.join(char for char in tissue if not char.isdigit())
            if "R" in tissue:
                tissue = "RL"   # Only using coarse grained annotations for right lung
            name = parts[1]
            return tissue, name
    return None, tip_name


def get_tissue_colors(tissues):
    """
    Generate distinct colors for each tissue type
    """
    # Manually set colors for specific tissue types
    manual_colors = {
        "LL": "black",
        "RL": "blue",
        "Liv": "green",
        "M": "red",
        "unknown": "lightgray"
    }
    
    # Additional colors for any other tissues that might be present
    fallback_colors = [
        "#FFEAA7", "#DDA0DD", "#98D8C8", "#F7DC6F", "#BB8FCE", 
        "#85C1E9", "#F8C471", "#82E0AA", "#F1948A", "#D7BDE2"
    ]
    
    tissue_colors = {}
    fallback_idx = 0
    
    for tissue in tissues:
        if tissue in manual_colors:
            tissue_colors[tissue] = manual_colors[tissue]
        else:
            # Use fallback colors for any other tissues
            tissue_colors[tissue] = fallback_colors[fallback_idx % len(fallback_colors)]
            fallback_idx += 1
    
    return tissue_colors


def make_tree_ultrametric(tree):
    """
    Make tree ultrametric by ensuring all tips are at the same distance from root
    """
    # Set all branch lengths to 1 for even spacing
    for node in tree.traverse():
        node.dist = 1
        
    # First, get the maximum distance from root to any tip
    max_distance = 0
    for leaf in tree.iter_leaves():
        distance = leaf.get_distance(tree)
        max_distance = max(max_distance, distance)
    
    # Now adjust all branch lengths so that all tips are at max_distance from root
    for leaf in tree.iter_leaves():
        current_distance = leaf.get_distance(tree)
        if current_distance < max_distance:
            # Add the difference to the leaf's branch length
            leaf.dist += (max_distance - current_distance)
    
    # Verify ultrametricity
    distances = [leaf.get_distance(tree) for leaf in tree.iter_leaves()]
    if len(set(distances)) > 1:
        print(f"Warning: Tree may not be perfectly ultrametric. Distance range: {min(distances)} - {max(distances)}")
    else:
        print(f"Tree is ultrametric. All tips at distance {distances[0]} from root")


def plot_circular_tree_by_tissue(newick_file, outfile):
    """
    Plot ultrametric tree in circular style with tips colored by tissue
    """
    # Read the tree
    tree = Tree(newick_file, format=1)
    
    # Make tree ultrametric for visualization
    make_tree_ultrametric(tree)
    
    # Extract tissues from tip names and create tissue mapping
    tissues = set()
    tissue_mapping = {}
    
    for node in tree.traverse():
        tissue, name = extract_tissue_from_name(node.name)
        if tissue:
            tissues.add(tissue)
            tissue_mapping[node] = tissue
            # Store the original name without tissue prefix
            node.original_name = name
        else:
            # If no tissue prefix found, use a default
            tissue_mapping[node] = "unknown"
            tissues.add("unknown")
            node.original_name = node.name
    
    # Convert to sorted list for consistent color assignment
    tissues = sorted(list(tissues))
    tissue_colors = get_tissue_colors(tissues)
    
    # Create tree style for circular layout
    ts = TreeStyle()
    ts.mode = "c"  # Circular mode
    ts.show_leaf_name = False  # Don't show tip names
    ts.show_branch_length = False  # Don't show branch lengths
    ts.show_scale = False  # Don't show scale
    ts.show_border = False  # No border
    
    # Add margins to prevent cutting off labels
    ts.margin_left = 10
    ts.margin_right = 10
    ts.margin_top = 10
    ts.margin_bottom = 10
    
    # Set up node styles
    for node in tree.traverse():
        nstyle = NodeStyle()

        # Color leaf nodes by tissue
        tissue = tissue_mapping.get(node, "unknown")
        color = tissue_colors[tissue]
        nstyle["fgcolor"] = color
        if node.is_leaf():
            nstyle["size"] = 8
        else:
            nstyle["size"] = 0
        nstyle["hz_line_color"] = color
        nstyle["vt_line_color"] = color
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_width"] = 1
        
        node.set_style(nstyle)
    
    # Add legend
    for tissue in tissues:
        color = tissue_colors[tissue]
        ts.legend.add_face(
            CircleFace(10, color, "circle"), 
            column=0
        )
        ts.legend.add_face(
            TextFace(f" {tissue}", fsize=10), 
            column=1
        )
    
    # Set environment for headless rendering
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    
    # Render the tree
    tree.render(outfile, tree_style=ts, dpi=600)
    print(f"Tree saved to {outfile}")
    print(f"Found {len(tissues)} tissue types: {', '.join(tissues)}")


def main():
    
    newick_file = sys.argv[1]
    outfile = sys.argv[2]

    # newick_file = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/mach2/5k/28/laml_trees_no_branch_lengths_no_origin.nwk"
    # outfile = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/mach2/5k/28/laml_trees_no_branch_lengths_no_origin.pdf"
    
    # Ensure output directory exists
    outdir = os.path.dirname(outfile)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
    
    plot_circular_tree_by_tissue(newick_file, outfile)


if __name__ == "__main__":
    main() 