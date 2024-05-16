#!/bin/bash

indel_matrix_path=$1
tissues_path=$2
template_xml=$3
REPLACE_TIME=$4
output_dir=$5

dir=$(dirname "$indel_matrix_path")
XML_FILE="$output_dir/joint_inference_beast.xml"
cp $template_xml $XML_FILE

# tissue labels for tips
REPLACE_TRAITSET=""
REPLACE_DATE_TRAITSET=""
while IFS= read -r line; do
    trait=$(echo $line | awk '{print $1 "=" $2 ","}')
    REPLACE_TRAITSET+="$trait"
    time=$(echo $line | awk '{print $1 "=REPLACE_TIME" ","}')
    REPLACE_DATE_TRAITSET+="$time"
done < "$tissues_path"

# total number of tissues
REPLACE_NUM_TISSUES=$(awk '{print $2}' "$tissues_path" | sort -u | wc -l)

# uniform tissue freqs
REPLACE_TISSUE_FREQUENCIES=$(echo "scale=10; 1 / $REPLACE_NUM_TISSUES" | bc)

# root tissue freqs
REPLACE_ROOT_TISSUE_FREQUENCIES="1"
for (( i=0; i < (( $REPLACE_NUM_TISSUES - 1 )); i++ )); do
REPLACE_ROOT_TISSUE_FREQUENCIES+=" 0"
done

# code map
tissues_string=$(awk '{print $2}' "$tissues_path" | sort -u)
tissues=($tissues_string)
primary_tissue="P"
non_primary_tissues=()
for item in "${tissues[@]}"; do
    if [[ $item != "$primary_tissue" ]]; then
        non_primary_tissues+=("$item")
    fi
done
sorted_np=($(for element in "${non_primary_tissues[@]}"; do echo "$element"; done | sort))
REPLACE_CODE_MAP="${primary_tissue}=0"
trailing_code_map=",? = 0"
for ((i=1; i<$REPLACE_NUM_TISSUES; i++)); do
    index=$((i-1))
    REPLACE_CODE_MAP+=",${sorted_np[index]}=$i"
    trailing_code_map+=" ${i}"
done
REPLACE_CODE_MAP+="${trailing_code_map}"

# num muts
unique_muts=$(tail -n +2 "$indel_matrix_path" | awk '{for (i=2; i<=NF; i++) print $i}' | sed 's/-1/0/g' | sort -u)
# ensure unedited state is included for edge case when all sites are edited and no dropouts occur since replacing dropouts as unedited here temporarily
if ! echo "$unique_muts" | grep -q '^0$'; then
    unique_muts+=" 0"
fi
REPLACE_NUM_MUTS=$(echo "$unique_muts" | wc -w)

REPLACE_NUM_EDITS=$((REPLACE_NUM_MUTS - 1))


REPLACE_SEQUENCES=$(python scripts/format_tidetree_sequences_sim_matrix.py $indel_matrix_path | tr -d '\n')
##### temp solution to not model site dropouts
REPLACE_SEQUENCES=$(printf '%s\n' "$REPLACE_SEQUENCES" | sed 's/-1/0/g')

# get the number of sites from the indel matrix
num_target_sites=$(tail -n 2 "$indel_matrix_path" | awk '{print NF-1; exit}' )
num_samples=$(tail -n +2 "$indel_matrix_path" | wc -l)
num_total_target_sites=$((num_target_sites * num_samples))


# calculated edit rates
all_muts=($(echo "$REPLACE_SEQUENCES" | grep -o "value='[0-9,]*'" | sed "s/value='//g" | sed "s/'//g" | sed "s/,/ /g" | tr '\n' ' '))
# Calculate value counts of all mutations
unset mut_counts
declare -A mut_counts
done_muts=()
for mut in "${all_muts[@]}"; do
    if [[ $mut == 0 ]]; then
        continue
    fi
    if [[ ! " ${done_muts[@]} " =~ " ${mut} " ]]; then
        done_muts+=("$mut")
    fi
    ((mut_counts[$mut]++))
done
# get the total count for the unedited state which will be set to 1.0 edit rate
ref_total=${mut_counts[${done_muts[0]}]}
editRates=()
# Print the mutation counts
for mut in ${done_muts[@]}; do
    editRate=("$(echo "scale=4; ${mut_counts[$mut]} / $ref_total" | bc)")
    formattedRate=$(printf "%.4f" "$editRate")
    editRates+=("$formattedRate")
done
REPLACE_EDIT_RATES=${editRates[@]}

# compute the rate needed to get the desired proportion_edited (computation is approximate for prior, not robustly calculated)
all_edits=$(echo "$REPLACE_SEQUENCES" | grep -o "value='[0-9,]*'" | sed "s/value='//g" | sed "s/'//g" | sed "s/,/ /g" | tr '\n' ' ' | sed 's/-1//g' | sed 's/0//g')
num_edits=$(echo "$all_muts" | wc -w)
proportion_edited=$(echo "scale=4; $num_edits / $num_total_target_sites" | bc)
# scaling since editRates are all 1.0 (again not the best approach but okay for now)
REPLACE_CLOCK_RATE=$(echo "scale=10; $proportion_edited / ($REPLACE_TIME)" | bc)

# format template xml file
REPLACE_SEQUENCES=$(printf '%s\n' "$REPLACE_SEQUENCES" | sed 's/[\/&]/\\&/g' | sed 's/\n//g' )
REPLACE_SEQUENCES="${REPLACE_SEQUENCES//\'/\"}"
REPLACE_TRAITSET=$(printf '%s\n' "$REPLACE_TRAITSET" | sed 's/[\/&]/\\&/g')

sed -i "s|REPLACE_SEQUENCES|$REPLACE_SEQUENCES|g" $XML_FILE
sed -i "s|REPLACE_TRAITSET|$REPLACE_TRAITSET|g" $XML_FILE
sed -i "s|REPLACE_DATE_TRAITSET|$REPLACE_DATE_TRAITSET|g" $XML_FILE
sed -i "s|REPLACE_NUM_TISSUES|$REPLACE_NUM_TISSUES|g" $XML_FILE
sed -i "s|REPLACE_TISSUE_FREQUENCIES|$REPLACE_TISSUE_FREQUENCIES|g" $XML_FILE
sed -i "s|REPLACE_CODE_MAP|$REPLACE_CODE_MAP|g" $XML_FILE
sed -i "s|REPLACE_ROOT_TISSUE_FREQUENCIES|$REPLACE_ROOT_TISSUE_FREQUENCIES|g" $XML_FILE
sed -i "s|REPLACE_NUM_MUTS|$REPLACE_NUM_MUTS|g" $XML_FILE
sed -i "s|REPLACE_NUM_EDITS|$REPLACE_NUM_EDITS|g" $XML_FILE
sed -i "s|REPLACE_EDIT_RATES|$REPLACE_EDIT_RATES|g" $XML_FILE
sed -i "s|REPLACE_TIME|$REPLACE_TIME|g" $XML_FILE
sed -i "s|REPLACE_CLOCK_RATE|$REPLACE_CLOCK_RATE|g" $XML_FILE

