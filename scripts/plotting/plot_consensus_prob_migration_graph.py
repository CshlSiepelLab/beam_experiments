#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import networkx as nx

DEFAULT_COLORS = ["#006400", "#FF0000", "#0000CD", "#FFA500", "#800080", "#808080", "#FFC0CB", "#ADD8E6", "#A52A2A", "#FFFF00"]*3


# inputs
graph_posterior_csv = sys.argv[1]
primary_tissue=sys.argv[2]
outdir = sys.argv[3]
consensus_probability_threshold = float(sys.argv[4])

# # testing
# graph_posterior_csv = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/metastabayes/MMUS1457/CP01/posterior_prob_graph.csv"
# primary_tissue="PRL"
# outdir = "."
# consensus_probability_threshold = 0.7


# obtain the probabilistic consensus migration graph
graph_dict = {}
with open(graph_posterior_csv, "r") as file:
    for line in file:
        key, value = line.strip().split(",")
        graph_dict[key] = float(value)

# find all tissues to set the node colors
all_tissues = sorted(list(set([value for node in graph_dict.keys() for value in node.split("_")[0:2]])))
custom_colors = DEFAULT_COLORS
custom_colors = {node: color for node, color in zip(all_tissues, custom_colors[0:len(all_tissues)]) if node != primary_tissue}
custom_colors[primary_tissue] = "black"

# plot the probability graph with edge thicknesses proportional to the probability
G = nx.MultiDiGraph()
for node in all_tissues:
    G.add_node(node, color=custom_colors[node], shape="box", fillcolor="white", penwidth=3.0, fontsize=32)
for edge, probability in graph_dict.items():
    source, target, num = edge.split('_')
    G.add_edge(source, target, color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"', penwidth=probability*3, fontsize=24)
dot = nx.nx_pydot.to_pydot(G)
dot.write_pdf(f"{outdir}/probability_migration_graph.pdf")

# plot the thresholded graph with all edges above consensus_probability_threshold probability
G = nx.MultiDiGraph()
for node in all_tissues:
    G.add_node(node, color=custom_colors[node], shape="box", fillcolor="white", penwidth=3.0, fontsize=32)
for edge, probability in graph_dict.items():
    if probability > consensus_probability_threshold:
        source, target, num = edge.split('_')
        G.add_edge(source, target, color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"', penwidth=3, fontsize=24)
dot = nx.nx_pydot.to_pydot(G)
dot.write_pdf(f"{outdir}/threshold_{consensus_probability_threshold}_migration_graph.pdf")

# plot the thresholded graph with all edges above consensus_probability_threshold probability but collapse directed multiedges into one with a number label of the original number of those edges
G = nx.MultiDiGraph()
for node in all_tissues:
    G.add_node(node, color=custom_colors[node], shape="box", fillcolor="white", penwidth=3.0, fontsize=32)
for edge, probability in graph_dict.items():
    if probability > consensus_probability_threshold:
        source, target, num = edge.split('_')
        if G.has_edge(source, target):
            G[source][target][0]['label'] = str(int(G[source][target][0]['label']) + 1)
        else:
            G.add_edge(source, target, color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"', penwidth=3, label="1", fontsize=24)

# Set the label to "" for any labels with just 1
for source, target, data in G.edges(data=True):
    if data.get('label') == '1':
        data['label'] = ''
        
dot = nx.nx_pydot.to_pydot(G)
dot.write_pdf(f"{outdir}/threshold_{consensus_probability_threshold}_migration_graph_collapsed_numbered.pdf")
