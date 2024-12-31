#!/bin/bash

# Paste normal BEAM setup here
##### BEAM setup for MMUS/CP combination #####

# calculate necessary inputs for the editing model
initialEditRates=$(cat /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1544/CP01/CP01_edit_rate_proportions.txt | tr '\n' ' ')
array=($initialEditRates)
numEditRates=${#array[@]}
numEditRatesPlusTwo=$(( numEditRates + 2 ))

# format tissue information
unique_tissues=$(cat /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1544/CP01/CP01_tip_tissues.csv | cut -d',' -f2 | sort | uniq)
numTissues=$(echo $unique_tissues | tr ' ' '\n' | wc -l)

# set the mcmc chain length based on the number of tips in the tree
# use a minimum of 5 million and a maximum of a user specified length
numTips=$(grep -c '^>' /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1544/CP01/CP01.fasta)
mcmclength=$(( numTips * 250000 ))
if (( mcmclength > 100000000 )); then
    mcmclength=100000000
elif (( mcmclength < 5000000 )); then
    mcmclength=5000000
fi

# mcmclength=200000000

# if numTissues is 1, then the model will not run, so we need to add a dummy tissue
if [[ $numTissues -eq 1 ]]; then
    # if the primary already exists then add a dummy
    if echo "$unique_tissues" | grep -w "PRL"; then
        unique_tissues="$unique_tissues dummy"
    # if not, then add the primary as the dummy
    else
        unique_tissues="PRL $unique_tissues"
    fi
    numTissues=$(( numTissues + 1 ))
fi

# force the known known primary tissue at the origin for root frequencies in the likelihood calculation
sorted_unique_tissues="PRL $(echo $unique_tissues | tr ' ' '\n' | grep -v '^PRL$' | tr '\n' ' ')"
if ! echo "$unique_tissues" | grep -qw "PRL"; then
    numTissues=$(( numTissues + 1 ))
fi
tissueRootFreqs=$(for i in $(seq 1 $numTissues); do if [[ "$i" = "1" ]]; then echo -n "1 "; else echo -n "0 "; fi; done)

numTissueRates=$(( numTissues * (numTissues - 1) ))
equalTissueRates=$(echo "scale=10; 1.0 / $numTissueRates" | bc -l)
equalTissueFreqs=$(echo "scale=10; 1.0 / $numTissues" | bc -l)

tissueCodeMap=""
trailingCodeMap="? = "
i=0
for tissue in $sorted_unique_tissues; do
    tissueCodeMap+="$tissue=$i,"
    trailingCodeMap+="$i "
    i=$((i+1))
done
tissueCodeMap+="${trailingCodeMap}"

# use LAML starting tree
# remove node names from newick and the origin node
newick=$(cat /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/laml/MMUS1544/CP01/laml_trees.nwk | sed 's/node[0-9]*:/:/g' | sed 's/^.//' | sed 's/..$//')


#####

mcmclength=1000000

# Run beam

outdir="test_beam"
beam_template="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/beam_mcmc_template_bd.xml"

cp $beam_template $outdir/beam.xml

# for some reason this cannot be taken in at the command line, likely because the = signs are problematic and not taken as a single string
sed -i "s/\$(tissueCodeMap)/$tissueCodeMap/g" $outdir/beam.xml


java -Xmx10g -jar /grid/siepel/home_norepl/staklins/beam/beam.jar \
    -threads 5 \
    -overwrite \
    -working \
    -D inName=CP01 \
    -D fileDir=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1544/CP01 \
    -D newick="$newick" \
    -D generations=224 \
    -D numEditRates="$numEditRates" \
    -D numEditRatesPlusTwo="$numEditRatesPlusTwo" \
    -D initialEditRates="$initialEditRates" \
    -D numTissueRates="$numTissueRates" \
    -D numTissues="$numTissues" \
    -D equalTissueFreqs="$equalTissueFreqs" \
    -D tissueRootFreqs="$tissueRootFreqs" \
    -D mcmclength="$mcmclength" \
    -D outname=test_beam \
    $outdir/beam.xml \
    > $outdir/beam.log 2>&1


# Call tidetree template with the same parameter inputs as required
# Update specific parameter inputs as needed

tidetree_xml_template="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/tidetree_template.xml"
outdir="test_tidetree"

cp $tidetree_xml_template $outdir/tidetree.xml

java -Xmx10g -jar /grid/siepel/home_norepl/staklins/beam/beam.jar \
    -threads 5 \
    -overwrite \
    -working \
    -D inName=CP01 \
    -D fileDir=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1544/CP01 \
    -D newick="$newick" \
    -D generations=224 \
    -D numEditRates="$numEditRates" \
    -D numEditRatesPlusTwo="$numEditRatesPlusTwo" \
    -D initialEditRates="$initialEditRates" \
    -D mcmclength="$mcmclength" \
    -D outname=test_tidetree \
    $outdir/tidetree.xml \
    > $outdir/tidetree.log 2>&1
