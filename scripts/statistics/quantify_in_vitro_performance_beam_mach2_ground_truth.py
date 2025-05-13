#!/usr/bin/env python3

import sys
import pandas as pd

inferred = sys.argv[1]
ground_truth = sys.argv[2]
outfile = sys.argv[3]

# Read the files into pandas DataFrames
inferred_df = pd.read_csv(inferred, names=["event", "prob"])
ground_truth_df = pd.read_csv(ground_truth)

# Setup source and target tissues for inferred and mach2 data
inferred_df["source"] = inferred_df["event"].str.split("_").str[0]
inferred_df["target"] = inferred_df["event"].str.split("_").str[1:-1].str.join(".")
inferred_df["count"] = inferred_df["event"].str.split("_").str[-1]
inferred_df = inferred_df.drop(columns=["event"])
inferred_df['prob'] = inferred_df['prob'].astype(float)
inferred_df['count'] = inferred_df['count'].astype(int)

# Since the ground truth multi-edge number is not really known in the same sense at the inference methods used, only use count 1 edges for inferred and mach2
ground_truth_df = ground_truth_df[ground_truth_df['count'] == 1]
inferred_df = inferred_df[inferred_df['count'] == 1]

# Subset the ground truth dataframe to only include the tissues possible in this cp
all_tissues = set(inferred_df['source'].unique()) | set(inferred_df['target'].unique())
ground_truth_df = ground_truth_df[ground_truth_df['source'].isin(all_tissues) & ground_truth_df['target'].isin(all_tissues)]
ground_truth_edge_set = set(zip(ground_truth_df['source'], ground_truth_df['target'], ground_truth_df['count']))

# Set probability thresholds to evaluate performance at
thresholds = [round(x * 0.01, 2) for x in range(0, 100)]

# Initialize DataFrames for metrics
metrics_inferred = pd.DataFrame(columns=['threshold', 'true_positives', 'false_positives', 'accuracy', 'precision', 'f1'])

# Check each edge in ground truth against inferred and mach2
for threshold in thresholds:
    # subset inferred and mach2 dataframe to prob >= threshold
    inferred_subset = inferred_df[inferred_df['prob'] >= threshold]
    
    inferred_edge_set = set(zip(inferred_subset['source'], inferred_subset['target'], inferred_subset['count']))
    
    # Calculate metrics for inferred
    tp_inferred = len(inferred_edge_set & ground_truth_edge_set)
    fp_inferred = len(inferred_edge_set - ground_truth_edge_set)

    accuracy = tp_inferred / len(ground_truth_edge_set)
    precision = tp_inferred / (tp_inferred + fp_inferred) if (tp_inferred + fp_inferred) > 0 else 0
    f1 = 2 * (precision * accuracy) / (precision + accuracy) if (precision + accuracy) > 0 else 0
    
    # Create a new row for this threshold
    new_row = {
        'threshold': threshold,
        'true_positives': tp_inferred,
        'false_positives': fp_inferred,
        'accuracy': accuracy,
        'precision': precision,
        'f1': f1
    }
    
    # Add the row to the DataFrame
    metrics_inferred = pd.concat([metrics_inferred, pd.DataFrame([new_row])], ignore_index=True)

# Save metrics to csv output files
metrics_inferred.to_csv(outfile, index=False)
