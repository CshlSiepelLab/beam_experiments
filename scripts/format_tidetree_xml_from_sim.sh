#!/bin/bash

# This script takes in files containing the information which needs to be plugged into the template xml file provided by tidetree repo to run on new data. It provides automation of the process of manually formatting tidetree input xml file.
# Designed with simulated data in mind but should technically work for any proper input files


if [[ $# -eq 0 ]] ; then
    echo "Usage: format_tidetree_xml_from_sim.sh --seqs <formatted sequences xml (str)> --total <total experiment time in hours(int)> --edit <edit time in hours as subset of total time (int)> --chain <length of mcmc chain for beast2, use 1000000000 for real runs (int)>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--seqs) seq_file="$2"; shift ;;
        -t|--total) total_time="$2"; shift ;;
        -e|--edit) edit_time="$2"; shift ;;
        -c|--chain) chain_length="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: format_tidetree_xml_from_sim.sh --seqs <formatted sequences xml (str)> --total <total experiment time in hours(int)> --edit <edit time in hours as subset of total time (int)> --chain <length of mcmc chain for beast2, use 1000000000 for real runs (int)>"; exit 1 ;;
    esac
    shift
done

REPLACE_SEQUENCES=""
states_array=()
cell_ids=()
while IFS= read -r line; do
    REPLACE_SEQUENCES+="$line "
    
    # Append values to array to define future variables on number of states and edit rates
    values=$(echo "$line" | grep -oE "value='([^']+)" | sed 's/,$//' | awk -F"'" '{print $2}')
    IFS=',' read -ra values_array <<< "$values"
    states_array+=(${values_array[@]})
    
    # Append cell id's to array to use later with defining tip dates
    name=$(echo "$line" | grep -oE "id='([^']+)" | awk -F"'" '{print $2}')
    cell_ids+=(${name})
    
done < "$seq_file"
max=$(IFS=$'\n'; echo "${states_array[*]}" | sort -nr | head -n1)

REPLACE_NUM_STATES=$((max + 1))

# Hack to calculate num sequences, but could also be user manual input if needed
NUM_CELLS=$(echo "$REPLACE_SEQUENCES" | grep -o '<' | wc -l)

# This is an integer number of hours
REPLACE_TOTAL_TIME=${total_time}

# This is an integer number of hours
REPLACE_EDIT_TIME=${edit_time}

# This is the MCMC chain length for beast. Tidetree default is 1000000000
REPLACE_CHAIN_LENGTH=$chain_length

# This is edit rates for one less than the number of states it seems (or number of columns in row 1 minue the first column which is default calculated)
# These priors are set to be equal rates for all scarring states, but can probably be better set in the future to represent expected rates or based on
# TideTree example had set bias for first rate to be much higher than the others, but it is not clear that there is a basis for this decision and the chain should run long enough to converge to the correct value regardless
equal_rate_value=$(echo "scale=2; 1 / ($REPLACE_NUM_STATES - 1)" | bc)
last_rate_value=$(echo "scale=2; 1 - $equal_rate_value * ($REPLACE_NUM_STATES - 2)" | bc)
REPLACE_EDIT_RATES=""
for ((i = 1; i < REPLACE_NUM_STATES; i++)); do
        if [ $i -eq $((REPLACE_NUM_STATES - 1)) ]; then
        REPLACE_EDIT_RATES+="0$last_rate_value"
    else
        REPLACE_EDIT_RATES+="0$equal_rate_value "
    fi
done

# this is the number of taxon/cells and total time for each cell
REPLACE_TIP_DATES=""
for ((iterator=0; iterator<NUM_CELLS; iterator++)); do
    REPLACE_TIP_DATES+="${cell_ids[iterator]}=$REPLACE_TOTAL_TIME"

    # Add a comma if it's not the last iteration
    if [ $((iterator+1)) -lt $NUM_CELLS ]; then
        REPLACE_TIP_DATES+=","
    fi
done


# This is 1 for first state (assumed to be unedited) and then 0 for the rest up to the number of possible states with space between each frequency
REPLACE_STARTING_FREQUENCIES="1"
for ((iterator=2; iterator<=REPLACE_NUM_STATES; iterator++)); do
    REPLACE_STARTING_FREQUENCIES+=" 0"
done

# Copy tidetree template to modify
XML_FILE=$(echo "$seq_file" | sed 's/\.xml/_formatted_for_tidetree.xml/')
cp inputs/tidetree_template.xml $XML_FILE

# Replace key words
REPLACE_SEQUENCES=$(printf '%s\n' "$REPLACE_SEQUENCES" | sed 's/[\/&]/\\&/g')
sed -i '' "s/REPLACE_SEQUENCES/$REPLACE_SEQUENCES/g" $XML_FILE
sed -i '' "s/REPLACE_NUM_STATES/$REPLACE_NUM_STATES/g" $XML_FILE
sed -i '' "s/REPLACE_TOTAL_TIME/$REPLACE_TOTAL_TIME/g" $XML_FILE
sed -i '' "s/REPLACE_EDIT_TIME/$REPLACE_EDIT_TIME/g" $XML_FILE
sed -i '' "s/REPLACE_CHAIN_LENGTH/$REPLACE_CHAIN_LENGTH/g" $XML_FILE
sed -i '' "s/REPLACE_EDIT_RATES/$REPLACE_EDIT_RATES/g" $XML_FILE
sed -i '' "s/REPLACE_TIP_DATES/$REPLACE_TIP_DATES/g" $XML_FILE
sed -i '' "s/REPLACE_STARTING_FREQUENCIES/$REPLACE_STARTING_FREQUENCIES/g" $XML_FILE


