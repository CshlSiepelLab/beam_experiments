#!/bin/bash

# This script takes in files containing the information which needs to be plugged into the template xml file provided by tidetree repo to run on new data. It provides automation of the process of manually formatting tidetree input xml file.
# Designed with simulated data in mind but should technically work for any proper input files

# This can be automated to read in lines from a file, but leaving manual for testing
# Need some whitespace between sequences, ideally would be newlines for formatting but space is used here for testing
#seq_file=$1
#total_time=$2
#edit_time=$3
#output_path=$4

seq_file="inputs/tidetree_seqs_example.xml"
total_time=54
edit_time=36
output_path="tidetree.xml"

REPLACE_SEQUENCES=""
states_array=()
while IFS= read -r line; do
    REPLACE_SEQUENCES+="$line "
    values=$(echo "$line" | grep -oE "value='([^']+)" | sed 's/,$//' | awk -F"'" '{print $2}')
    IFS=',' read -ra values_array <<< "$values"
    states_array+=(${values_array[@]})
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
REPLACE_CHAIN_LENGTH=1000000

# This is edit rates for one less than the number of states it seems
REPLACE_EDIT_RATES="0.9 0.1"

# this is the number of taxon/cells and total time for each cell
REPLACE_TIP_DATES=""
for ((iterator=1; iterator<=NUM_CELLS; iterator++)); do
    REPLACE_TIP_DATES+="$iterator=$REPLACE_TOTAL_TIME"

    # Add a comma if it's not the last iteration
    if [ $iterator -lt $NUM_CELLS ]; then
        REPLACE_TIP_DATES+=","
    fi
done


# This is 1 for first state (assumed to be unedited) and then 0 for the rest up to the number of possible states with space between each frequency
REPLACE_STARTING_FREQUENCIES="1"
for ((iterator=2; iterator<=REPLACE_NUM_STATES; iterator++)); do
    REPLACE_STARTING_FREQUENCIES+=" 0"
done

# User specified output dir set as is for testing purposes, but should be modified to be input
OUTDIR=""

# Copy tidetree template to modify
XML_FILE=${output_path}
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


