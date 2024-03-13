#!/bin/bash

# This script takes in files containing the information which needs to be plugged into the template xml file for FixedTreeAnalysis.

# seqfile=$1
# taxafile=$2
# traitfile=$3
# newickfile=$4
# xml_template=$5
# primary_tissue=$6
# chainlength=$7
# model=$8    ### can be "asym", "sym", "oneRate", "threeRates"

seqfile="results/fixed_offset_sym_asym_compare_machina_sims_fixedtreeanalysis_2_21_24/fixed_offset_symmmetricalfalse_machina_m8_sims_compare_beast_machina_fixedtreeanalysis_default_2_21_24/machina_m8_sim_data/seed172/T_seed172_unlabeled_true_tree_sequences_formatted_for_xml.txt"
taxafile="results/fixed_offset_sym_asym_compare_machina_sims_fixedtreeanalysis_2_21_24/fixed_offset_symmmetricalfalse_machina_m8_sims_compare_beast_machina_fixedtreeanalysis_default_2_21_24/machina_m8_sim_data/seed172/T_seed172_unlabeled_true_tree_taxonset_formatted_for_xml.txt"
traitfile="results/fixed_offset_sym_asym_compare_machina_sims_fixedtreeanalysis_2_21_24/fixed_offset_symmmetricalfalse_machina_m8_sims_compare_beast_machina_fixedtreeanalysis_default_2_21_24/machina_m8_sim_data/seed172/T_seed172_unlabeled_true_tree_traitset_formatted_for_xml.txt"
newickfile="results/fixed_offset_sym_asym_compare_machina_sims_fixedtreeanalysis_2_21_24/fixed_offset_symmmetricalfalse_machina_m8_sims_compare_beast_machina_fixedtreeanalysis_default_2_21_24/machina_m8_sim_data/seed172/T_seed172_unlabeled_true_tree_newick_formatted_for_xml.txt"
xml_template="inputs/no_bsvss_template_xml_fixedtreeanalysis_machina_sim_universal.xml"
primary_tissue="P"
model="oneRate"

REPLACE_CHAINLENGTH="${chainlength}"

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

traits=()
while IFS= read -r line; do
    REPLACE_TRAITSET+="$line "
    trait=$(echo $line | awk -F'=' '{print $2}' | awk -F',' '{print $1}')
    if [[ ! " ${traits[@]} " =~ " $trait " ]]; then
        traits+=("$trait")
    fi
done < "$traitfile"
REPLACE_NUM_TISSUES=${#traits[@]}
REPLACE_OFFSET=$(( $REPLACE_NUM_TISSUES - 1 ))
REPLACE_TISSUE_FREQS=$(echo "scale=10; 1 / $REPLACE_NUM_TISSUES" | bc)

if [ "$model" = "sym" ]; then
    REPLACE_NUM_RATES=$(((REPLACE_NUM_TISSUES * (REPLACE_NUM_TISSUES - 1)) / 2))
    REPLACE_SUBSTITUTION_MODEL="<substModel id=\"svs.s:tissue\" spec=\"beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel\" rateIndicator=\"true\" rates=\"@relativeGeoRates.s:tissue\" symmetric=\"true\">"
    REPLACE_GEO_LOGGER="<log id=\"geoSubstModelLogger.s:tissue\" spec=\"beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModelLogger\" dataType=\"@traitDataType.tissue\" model=\"@svs.s:tissue\"/>"
elif [ "$model" = "asym" ]; then
    REPLACE_NUM_RATES=$((REPLACE_NUM_TISSUES * (REPLACE_NUM_TISSUES - 1)))
    REPLACE_SUBSTITUTION_MODEL="<substModel id=\"svs.s:tissue\" spec=\"beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel\" rateIndicator=\"true\" rates=\"@relativeGeoRates.s:tissue\" symmetric=\"false\">"
    REPLACE_GEO_LOGGER="<log id=\"geoSubstModelLogger.s:tissue\" spec=\"beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModelLogger\" dataType=\"@traitDataType.tissue\" model=\"@svs.s:tissue\"/>"
elif [ "$model" = "oneRate" ]; then
    REPLACE_NUM_RATES=1
    REPLACE_SUBSTITUTION_MODEL="<substModel id=\"svs.s:tissue\" spec=\"metastabayes.substitutionmodel.OneRateAllTissues\" rates=\"@relativeGeoRates.s:tissue\" rateIndicator=\"true\" symmetric=\"true\">"
    REPLACE_GEO_LOGGER="<log id=\"geoSubstModelLogger.s:tissue\" spec=\"metastabayes.substitutionmodel.OneRateAllTissuesLogger\" dataType=\"@traitDataType.tissue\" model=\"@svs.s:tissue\"/>"
elif [ "$model" = "threeRate" ]; then
    REPLACE_NUM_RATES=3
    REPLACE_SUBSTITUTION_MODEL="<substModel id=\"svs.s:tissue\" spec=\"metastabayes.substitutionmodel.ThreeRatesForSeedingRoutes\" rates=\"@relativeGeoRates.s:tissue\" rateIndicator=\"true\" symmetric=\"false\">"
    REPLACE_GEO_LOGGER="<log id=\"geoSubstModelLogger.s:tissue\" spec=\"metastabayes.substitutionmodel.ThreeRatesForSeedingRoutesLogger\" dataType=\"@traitDataType.tissue\" model=\"@svs.s:tissue\"/>"
else
    echo "Model input flag does not match allowable options. Please use asym, sym, oneRate, or threeRates."
    exit
fi

non_primary_traits=()
for item in "${traits[@]}"; do
    if [[ $item != "$primary_tissue" ]]; then
        non_primary_traits+=("$item")
    fi
done

sorted_np=($(for element in "${non_primary_traits[@]}"; do echo "$element"; done | sort))

REPLACE_CODE_MAP="${primary_tissue}=0"
trailing_code_map=",? = 0"
REPLACE_ROOT_FREQUENCIES="1"
for ((i=1; i<$REPLACE_NUM_TISSUES; i++)); do
    REPLACE_ROOT_FREQUENCIES+=" 0"
    index=$((i-1))
    REPLACE_CODE_MAP+=",${sorted_np[index]}=$i"
    trailing_code_map+=" ${i}"
done

REPLACE_CODE_MAP+="${trailing_code_map}"


while IFS= read -r line; do
    REPLACE_NEWICK+="$line"
done < "$newickfile"
# Copy xml template to modify
XML_FILE=$(echo "$seqfile" | sed "s/\_sequences_formatted_for_xml.txt/_final_input_xml_${model}.xml/")
cp $xml_template $XML_FILE


REPLACE_SEQUENCES=$(printf '%s\n' "$REPLACE_SEQUENCES" | sed 's/[\/&]/\\&/g')
REPLACE_TAXONSET=$(printf '%s\n' "$REPLACE_TAXONSET" | sed 's/[\/&]/\\&/g')
REPLACE_TRAITSET=$(printf '%s\n' "$REPLACE_TRAITSET" | sed 's/[\/&]/\\&/g')
REPLACE_NEWICK=$(printf '%s\n' "$REPLACE_NEWICK" | sed 's/[\/&]/\\&/g')

# Replace single quotes with double quotes for sequences and taxon
REPLACE_SEQUENCES="${REPLACE_SEQUENCES//\'/\"}"
REPLACE_TAXONSET="${REPLACE_TAXONSET//\'/\"}"

sed -i "s|REPLACE_SEQUENCES|$REPLACE_SEQUENCES|g" $XML_FILE
sed -i "s|REPLACE_SUBSTITUTION_MODEL|$REPLACE_SUBSTITUTION_MODEL|g" $XML_FILE
sed -i "s|REPLACE_GEO_LOGGER|$REPLACE_GEO_LOGGER|g" $XML_FILE
sed -i "s|REPLACE_TAXONSET|$REPLACE_TAXONSET|g" $XML_FILE
sed -i "s|REPLACE_TRAITSET|$REPLACE_TRAITSET|g" $XML_FILE
sed -i "s|REPLACE_NEWICK|$REPLACE_NEWICK|g" $XML_FILE
sed -i "s|REPLACE_NUM_TISSUES|$REPLACE_NUM_TISSUES|g" $XML_FILE
sed -i "s|REPLACE_TISSUE_FREQS|$REPLACE_TISSUE_FREQS|g" $XML_FILE
sed -i "s|REPLACE_NUM_RATES|$REPLACE_NUM_RATES|g" $XML_FILE
sed -i "s|REPLACE_CODE_MAP|$REPLACE_CODE_MAP|g" $XML_FILE
sed -i "s|REPLACE_ROOT_FREQUENCIES|$REPLACE_ROOT_FREQUENCIES|g" $XML_FILE
sed -i "s|REPLACE_CHAINLENGTH|$REPLACE_CHAINLENGTH|g" $XML_FILE
sed -i "s|REPLACE_OFFSET|$REPLACE_OFFSET|g" $XML_FILE


# add metastabayes package to namespace if necessary based on the input substitution model selection
if [ "$model" = "oneRate" ] || [ "$model" = "threeRates" ]; then
    new_namespace=":metastabayes"
    sed -i "1s/namespace=\"\(.*\)\"/namespace=\"\1$new_namespace\"/" "$XML_FILE"
fi

