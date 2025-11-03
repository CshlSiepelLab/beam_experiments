import sys
import ast
import math


metient_file = sys.argv[1]
outfile = sys.argv[2]

metient_prob_graph = {}
num_solutions = None

# Get all input graphs
with open(metient_file, "r") as file:
    graph_lines = file.readlines()[1:]  # skip header line
    num_solutions = len(graph_lines)
    all_graphs = []
    for line in graph_lines:
        loss, graph = line.strip().split("\t")
        all_graphs.append((float(loss), ast.literal_eval(graph)))

# Decide whether the consensus is equally weighted or weighted by loss
prob = None
metient_equal_prob = False  # Use if each sampled graph should have equal probability
if metient_equal_prob:
    prob = 1.0 / num_solutions
else:
    max_loss = max(solution[0] for solution in all_graphs)
    min_loss = min(solution[0] for solution in all_graphs)
    min_max_denominator = max_loss - min_loss if max_loss != min_loss else 1.0  # Prevent division by zero when all losses are equal
    all_graphs = [((loss - min_loss)/min_max_denominator, counts) for loss, counts in all_graphs]  # Min-max scale losses
    temp = 0.5  # Temperature parameter for softmax (fixed to 0.5 from Metient authors suggestion)
    prob_denominator = sum([math.exp(-loss/temp) for loss, counts in all_graphs])
    all_graphs = [(math.exp(-loss/temp)/prob_denominator, counts) for loss, counts in all_graphs]   # Convert to probabilities by temperature-scaled softmax
    
# Build the consensus graph
for solution_num, (loss, metient_counts_input) in enumerate(all_graphs):
    if not prob:
        prob = loss
    for source_tissue, targets_dict in metient_counts_input.items():
        for target_tissue, edge_count in targets_dict.items():
            if edge_count > 0:
                for n in range(1, int(edge_count) + 1):
                    migration = f"{source_tissue}_{target_tissue}_{n}"
                    if migration not in metient_prob_graph:
                        metient_prob_graph[migration] = prob
                    else:
                        metient_prob_graph[migration] += prob

# Output the consensus graph
with open(outfile, "w") as f:
    for migration, probability in metient_prob_graph.items():
        f.write(f"{migration},{probability}\n")
