import csv
import sys

input_csv = sys.argv[1]
output_dot = sys.argv[2]

# Assign tissue colors in order of appearance
color_map = {}
next_color = 1

def get_color(node):
    global next_color
    if node not in color_map:
        color_map[node] = next_color
        next_color += 1
    return color_map[node]

with open(input_csv) as f, open(output_dot, "w") as out:
    reader = csv.DictReader(f)

    out.write("digraph {\n")
    out.write("\tnode [colorscheme=set19 penwidth=3 shape=box]\n")
    out.write("\tedge [colorscheme=set19 penwidth=3]\n")

    # First pass: gather nodes
    edges = []
    for row in reader:
        source, target = row["source_target"].split("_")
        num_edges = int(row["num_edges"])
        edges.append((source, target, num_edges))
        get_color(source)
        get_color(target)

    # Write node lines
    for node, color in color_map.items():
        out.write(f'\t{node} [color={color}]\n')

    # Write edge lines
    for source, target, num_edges in edges:
        src_color = color_map[source]
        tgt_color = color_map[target]
        mix = f'"{src_color};0.5:{tgt_color}"'
        for _ in range(num_edges):
            out.write(f'\t{source} -> {target} [color={mix}]\n')

    out.write("}\n")
