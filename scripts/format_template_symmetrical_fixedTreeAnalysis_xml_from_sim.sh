#!/bin/bash

# This script takes in files containing the information which needs to be plugged into the template xml file for symmetrical rate matrix FixedTreeAnalysis.

seqfile=$1
taxafile=$2
traitfile=$3
newickfile=$4
xml_template=$5

REPLACE_SEQUENCES=""
REPLACE_NEWICK=""
REPLACE_TAXONSET=""
REPLACE_TRAITSET=""

while IFS= read -r line; do
    REPLACE_SEQUENCES+="$line "
done < "$seqfile"

while IFS= read -r line; do
    REPLACE_TAXONSET+="$line "
done < "$taxafile"

while IFS= read -r line; do
    REPLACE_TRAITSET+="$line "
done < "$traitfile"

while IFS= read -r line; do
    REPLACE_NEWICK+="$line"
done < "$newickfile"

# Copy xml template to modify
XML_FILE=$(echo "$seqfile" | sed 's/\_sequences_formatted_for_xml.txt/_final_input_xml.xml/')
cp $xml_template $XML_FILE


REPLACE_SEQUENCES=$(printf '%s\n' "$REPLACE_SEQUENCES" | sed 's/[\/&]/\\&/g')
REPLACE_TAXONSET=$(printf '%s\n' "$REPLACE_TAXONSET" | sed 's/[\/&]/\\&/g')
REPLACE_TRAITSET=$(printf '%s\n' "$REPLACE_TRAITSET" | sed 's/[\/&]/\\&/g')
REPLACE_NEWICK=$(printf '%s\n' "$REPLACE_NEWICK" | sed 's/[\/&]/\\&/g')

# Replace single quotes with double quotes for sequences and taxon
REPLACE_SEQUENCES="${REPLACE_SEQUENCES//\'/\"}"
REPLACE_TAXONSET="${REPLACE_TAXONSET//\'/\"}"

# Replace key words; Requires '' in sed command on mac, but not on linux
# sed -i '' "s|REPLACE_SEQUENCES|$REPLACE_SEQUENCES|g" $XML_FILE
# sed -i '' "s|REPLACE_TAXONSET|$REPLACE_TAXONSET|g" $XML_FILE
# sed -i '' "s|REPLACE_TRAITSET|$REPLACE_TRAITSET|g" $XML_FILE
# sed -i '' "s|REPLACE_NEWICK|$REPLACE_NEWICK|g" $XML_FILE

sed -i "s|REPLACE_SEQUENCES|$REPLACE_SEQUENCES|g" $XML_FILE
sed -i "s|REPLACE_TAXONSET|$REPLACE_TAXONSET|g" $XML_FILE
sed -i "s|REPLACE_TRAITSET|$REPLACE_TRAITSET|g" $XML_FILE
sed -i "s|REPLACE_NEWICK|$REPLACE_NEWICK|g" $XML_FILE
