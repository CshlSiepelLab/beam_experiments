#!/bin/bash

indel_matrix_path=$1
tissues_path=$2
template_xml=$3

dir=$(dirname "$indel_matrix_path")
XML_FILE="$dir/joint_inference_beast.xml"
cp $template_xml $XML_FILE

# tissue labels for tips
REPLACE_TRAITSET=""
REPLACE_DATE_TRAITSET=""
while IFS= read -r line; do
    trait=$(echo $line | awk '{print $1 "=" $2 ","}')
    REPLACE_TRAITSET+="$trait"
    time=$(echo $line | awk '{print $1 "=54" ","}')
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
REPLACE_NUM_MUTS=$(tail -n +2 "$indel_matrix_path" | awk '{for (i=2; i<=NF; i++) print $i}' | sort -u | wc -l)

# edit root frequencies
REPLACE_EDIT_ROOT_FREQUENCIES="1"
for (( i=0; i < (( $REPLACE_NUM_MUTS - 1 )); i++ )); do
REPLACE_EDIT_ROOT_FREQUENCIES+=" 0"
done

# sequences formatted for indel matrix
sequences=($(tail -n +2 "$indel_matrix_path" | cut -f1- |tr '\t' ','))
i=1
REPLACE_SEQUENCES=""
for seq in ${sequences[@]}; do
name=$(echo "$seq" | cut -d',' -f1)
muts=$(echo "$seq" | cut -d',' -f2-)
REPLACE_SEQUENCES+="<sequence id='${name}' spec='Sequence' taxon='${name}' value='${muts},'/> "
i=$(( i + 1 ))
done
##### temp solution to not model site dropouts
REPLACE_SEQUENCES=$(printf '%s\n' "$REPLACE_SEQUENCES" | sed 's/-1/0/g')

# edit rates
equal_rate_value=$(echo "scale=4; 1 / ($REPLACE_NUM_MUTS - 1)" | bc)
last_rate_value=$(echo "scale=4; 1 - $equal_rate_value * ($REPLACE_NUM_MUTS - 2)" | bc)
REPLACE_EDIT_RATES=""
for ((i = 1; i < REPLACE_NUM_MUTS; i++)); do
        if [ $i -eq $((REPLACE_NUM_MUTS - 1)) ]; then
        REPLACE_EDIT_RATES+="0$last_rate_value"
    else
        REPLACE_EDIT_RATES+="0$equal_rate_value "
    fi
done

# format template xml file
REPLACE_SEQUENCES=$(printf '%s\n' "$REPLACE_SEQUENCES" | sed 's/[\/&]/\\&/g')
REPLACE_TRAITSET=$(printf '%s\n' "$REPLACE_TRAITSET" | sed 's/[\/&]/\\&/g')
REPLACE_SEQUENCES="${REPLACE_SEQUENCES//\'/\"}"

sed -i "s|REPLACE_SEQUENCES|$REPLACE_SEQUENCES|g" $XML_FILE
sed -i "s|REPLACE_TRAITSET|$REPLACE_TRAITSET|g" $XML_FILE
sed -i "s|REPLACE_DATE_TRAITSET|$REPLACE_DATE_TRAITSET|g" $XML_FILE
sed -i "s|REPLACE_NUM_TISSUES|$REPLACE_NUM_TISSUES|g" $XML_FILE
sed -i "s|REPLACE_TISSUE_FREQUENCIES|$REPLACE_TISSUE_FREQUENCIES|g" $XML_FILE
sed -i "s|REPLACE_CODE_MAP|$REPLACE_CODE_MAP|g" $XML_FILE
sed -i "s|REPLACE_ROOT_TISSUE_FREQUENCIES|$REPLACE_ROOT_TISSUE_FREQUENCIES|g" $XML_FILE
sed -i "s|REPLACE_NUM_MUTS|$REPLACE_NUM_MUTS|g" $XML_FILE
sed -i "s|REPLACE_EDIT_ROOT_FREQUENCIES|$REPLACE_EDIT_ROOT_FREQUENCIES|g" $XML_FILE
sed -i "s|REPLACE_EDIT_RATES|$REPLACE_EDIT_RATES|g" $XML_FILE

