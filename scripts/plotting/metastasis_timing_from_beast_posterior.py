import pandas as pd
import seaborn as sns
from datetime import datetime
from collections import defaultdict
import dendropy

import matplotlib.pyplot as plt

def parse_tree(tree):
    migration_events = []
    for node in tree.postorder_node_iter():
        if node.edge_length is not None and node.annotations.get_value("location"):
            parent_location = node.parent_node.annotations.get_value("location") if node.parent_node else None
            current_location = node.annotations.get_value("location")
            if parent_location and parent_location != current_location:
                migration_events.append({
                    'start_time': node.parent_node.age,
                    'end_time': node.age,
                    'from': parent_location,
                    'to': current_location
                })
    return migration_events

def read_trees(file_path):
    trees = dendropy.TreeList.get(path=file_path, schema="nexus")
    return trees

def extract_migration_events(trees):
    all_events = []
    for tree in trees:
        events = parse_tree(tree)
        all_events.extend(events)
    return all_events


posterior_file = "/path/to/your/beast2_posterior_trees.nexus"


trees = read_trees(posterior_file)

migration_events = extract_migration_events(trees)

df = pd.DataFrame(migration_events)
df['start_time'] = pd.to_datetime(df['start_time'], unit='s')
df['end_time'] = pd.to_datetime(df['end_time'], unit='s')

plt.figure(figsize=(12, 8))
sns.set(style="whitegrid")

for i, event in df.iterrows():
    plt.plot([event['start_time'], event['end_time']], [event['from'], event['to']], marker='o')

plt.xlabel('Time')
plt.ylabel('Tissue Location')
plt.title('Migration Events Over Time')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()