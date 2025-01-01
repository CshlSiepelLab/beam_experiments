#!/bin/bash

# Paste normal BEAM setup here
##### BEAM setup for MMUS/CP combination #####

# calculate necessary inputs for the editing model
initialEditRates=$(cat /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_12_31_24_uniform_50cells_50sites_data_7_24_24/beam/mS_854/mS_854_edit_rate_proportions.txt)
array=($initialEditRates)
numEditRates=${#array[@]}
numEditRatesPlusTwo=$(( numEditRates + 2 ))

# format tissue information
unique_tissues=$(cat /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_12_31_24_uniform_50cells_50sites_data_7_24_24/beam/mS_854/mS_854_tip_tissues.csv | cut -d',' -f2 | sort | uniq)
numTissues=$(echo $unique_tissues | tr ' ' '\n' | wc -l)

# force the known known primary tissue at the origin for root frequencies in the likelihood calculation
sorted_unique_tissues="P $(echo $unique_tissues | tr ' ' '\n' | grep -v '^P$' | tr '\n' ' ')"
if ! echo "$unique_tissues" | grep -qw "P"; then
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

#####

mcmclength=15000000

# Run beam

# outdir="test_beam"
# mkdir -p $outdir
# beam_template="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/beam_mcmc_template_bd.xml"

# cp $beam_template $outdir/beam.xml

# # for some reason this cannot be taken in at the command line, likely because the = signs are problematic and not taken as a single string
# sed -i "s/\$(tissueCodeMap)/$tissueCodeMap/g" $outdir/beam.xml


# java -Xmx10g -jar /grid/siepel/home_norepl/staklins/beam/beam.jar \
#     -threads 5 \
#     -overwrite \
#     -working \
#     -D inName=CP01 \
#     -D fileDir=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1544/CP01 \
#     -D newick="$newick" \
#     -D generations=224 \
#     -D numEditRates="$numEditRates" \
#     -D numEditRatesPlusTwo="$numEditRatesPlusTwo" \
#     -D initialEditRates="$initialEditRates" \
#     -D numTissueRates="$numTissueRates" \
#     -D numTissues="$numTissues" \
#     -D equalTissueFreqs="$equalTissueFreqs" \
#     -D tissueRootFreqs="$tissueRootFreqs" \
#     -D mcmclength="$mcmclength" \
#     -D outname=test_beam \
#     $outdir/beam.xml \
#     > $outdir/beam.log 2>&1


# Call tidetree template with the same parameter inputs as required
# Update specific parameter inputs as needed

tidetree_xml_template="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/tidetree_template.xml"
outdir="test_tidetree_mS_854"
mkdir -p $outdir

cp $tidetree_xml_template $outdir/tidetree.xml

java -Xmx10g -jar /grid/siepel/home_norepl/staklins/beam/beam.jar \
    -threads 5 \
    -overwrite \
    -working \
    -D inName=mS_854 \
    -D fileDir=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_12_31_24_uniform_50cells_50sites_data_7_24_24/beam/mS_854 \
    -D newick="$newick" \
    -D generations=250 \
    -D numEditRates="$numEditRates" \
    -D numEditRatesPlusTwo="$numEditRatesPlusTwo" \
    -D initialEditRates="$initialEditRates" \
    -D mcmclength="$mcmclength" \
    -D outname=test_tidetree \
    $outdir/tidetree.xml \
    > $outdir/tidetree.log 2>&1
