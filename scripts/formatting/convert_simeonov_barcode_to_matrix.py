
import os
import sys
import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]   # should have _matrix.csv suffix, else change below

os.makedirs(os.path.dirname(output_file), exist_ok=True)

output_tissues_file = output_file.replace("_matrix.csv", "_tissues.txt")
output_mut_mapping_file = output_file.replace("_matrix.csv", "_mut_mapping.txt")

df = pd.read_csv(input_file, sep='\t', index_col=0)

df.index = ["barcode" + str(idx) for idx in df.index]

# Use cells per barcode to determine tissues
tissues = {}
for idx, row in df.iterrows():
    cells = [cell.strip() for cell in row['cells'].split(',')]
    tissue_set = set(cell.split('_')[0] for cell in cells)
    tissues[idx] = tissue_set
    
df = df.drop(columns=['cells'])

# Split barcode into sites
hmid_split = df['hmid'].str.split('-', expand=True)
num_cols = hmid_split.shape[1]
hmid_split.columns = [f"r{i+1}" for i in range(num_cols)]
df = pd.concat([df.drop(columns=['hmid']), hmid_split], axis=1)

# Replace UKNNOWN with -1 and NONE with 0
df = df.replace("UNKNOWN", "-1")
df = df.replace("NONE", "0")

# Now replace mut string with integers successively while also assigning multiple site mut str to the first site with other sites as -1
mut_mapping = {}
next_mut_id = 1
for idx, row in df.iterrows():
    prev_val= None
    
    for col in df.columns:
        val = row[col]
        # Skip unedited or missing data sites
        if val == "-1" or val == "0":
            prev_val = None
        
        # if val == prev_val, set to -1
        elif val == prev_val:
            df.at[idx, col] = -1
        elif val in mut_mapping:
            df.at[idx, col] = mut_mapping[val]
            prev_val = val
        else:
            mut_mapping[val] = next_mut_id
            df.at[idx, col] = next_mut_id
            next_mut_id += 1
            prev_val = val

# Write all to output csv files
with open(output_mut_mapping_file, "w") as f:
    f.write("successive_char_int,mut_str\n")
    for mut_str, mut_id in mut_mapping.items():
        f.write(f"{mut_id},{mut_str}\n")
        
with open(output_tissues_file, "w") as f:
    f.write("group_name,tissues\n")
    for barcode, tissue_set in tissues.items():
        f.write(f"{barcode},{';'.join(sorted(tissue_set))}\n")

df.to_csv(output_file, index=True, header=True)
                