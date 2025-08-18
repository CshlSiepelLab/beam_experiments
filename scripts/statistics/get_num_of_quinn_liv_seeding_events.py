
import pandas as pd

infile="/grid/siepel/home/staklins/stored_results/beam/latest_results/general_graph_stats_for_beam_paper_from_latest_runs_8_2_25/beam_all_results_8_18_25.csv"

data = pd.read_csv(infile)

data[['source', 'target', 'edgenum']] = data['source_target_edgenum'].str.split('_', expand=True)

cps_with_liv = set()
for cp in data['cp'].unique():
    cp_data = data[data['cp'] == cp]
    tissues = set(cp_data['source']).union(set(cp_data['target']))
    if "Liv" in tissues:
        cps_with_liv.add(cp)

data = data[data['cp'].isin(cps_with_liv)]

thresholds = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

results = {}
for threshold in thresholds:
    results[threshold] = {'LL': 0, 'RL': 0, 'M': 0}

for threshold in thresholds:
    threshold_data = data[data['probability'] > threshold]
    for cp in threshold_data['cp'].unique():
        cp_data = threshold_data[threshold_data['cp'] == cp]
        for tissue in ['LL', 'RL', 'M']:
            if ((cp_data['source'] == tissue) & (cp_data['target'] == 'Liv')).any():
                results[threshold][tissue] += 1


    