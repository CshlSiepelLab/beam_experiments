
countscript=/grid/siepel/home/staklins/projects/crispr_barcode/bayesian_phylogenetic_metastasis/scripts/formatting/count_phylogenetic_information_per_site.py

# # Variable rates data
# inDir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/data/variable_migration_and_mutation_rates_8_19_24/"
# outFile="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_2_25_25_data_from_8_19_24/phylogenetic_information_per_site_variable_rates.csv"

# echo "name,average_informative_characters_per_site" > $outFile

# matrixFiles=$(find $inDir -type f -name "*indel_character_matrix.tsv")

# for matrixFile in $matrixFiles; do
#     name=$(basename $(dirname $matrixFile))
#     val=$(python $countscript $matrixFile | grep "Average informative characters per site:" | cut -d " " -f 6)
#     echo "$name,$val" >> $outFile
# done

# # Quinn data
# inDir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/successive_raw_data/"
# outFile="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/phylogenetic_information_per_site_quinn.csv"

# echo "name,average_informative_characters_per_site" > $outFile

# matrixFiles=$(find $inDir -type f -name "*original_character_matrix.tsv")

# for matrixFile in $matrixFiles; do
#     name=$(basename $matrixFile | cut -d "_" -f 1)
#     val=$(python $countscript $matrixFile | grep "Average informative characters per site:" | cut -d " " -f 6)
#     echo "$name,$val" >> $outFile
# done

# Serio data
inDir="/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/beam"
outFile="/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/phylogenetic_information_per_site_serio2.csv"

echo "name,average_informative_characters_per_site" > $outFile

matrixFiles=$(find $inDir -type f -name "temp_matrix.tsv")

for matrixFile in $matrixFiles; do
    name=$(echo $matrixFile | rev | cut -d "/" -f2-3 | rev | sed 's/\//_/g')
    val=$(python $countscript $matrixFile | grep "Average informative characters per site:" | cut -d " " -f 6)
    echo "$name,$val" >> $outFile
done

# Remove serio CP00 rows from the file after running the above script
sed -i '/CP00/d' $outFile
