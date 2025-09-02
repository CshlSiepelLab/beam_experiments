
countscript=/grid/siepel/home/staklins/projects/crispr_barcode/bayesian_phylogenetic_metastasis/scripts/formatting/count_phylogenetic_information.py


# Variable rates data
inDir="/grid/siepel/home/staklins/stored_data/crispr_barcode_related_data/variable_migration_and_mutation_rates_8_19_24/"
outFile="/grid/siepel/home/staklins/stored_results/beam/latest_results/variable_migration_and_mutation_rates_data_8_19_24/phylogenetic_information_to_cell_ratio_variable_rates.csv"

echo "name,informative_muts_to_cell_ratio" > $outFile

matrixFiles=$(find $inDir -type f -name "*indel_character_matrix.tsv")

for matrixFile in $matrixFiles; do
    name=$(basename $(dirname $matrixFile))
    val=$(python $countscript $matrixFile | grep "Result:" | cut -d " " -f 2)
    echo "$name,$val" >> $outFile
done


# Quinn data
inDir="/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/successive_raw_data/"
outFile="/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/phylogenetic_information_to_cell_ratio_quinn.csv"

echo "name,informative_muts_to_cell_ratio" > $outFile

matrixFiles=$(find $inDir -type f -name "*original_character_matrix.tsv")

for matrixFile in $matrixFiles; do
    name=$(basename $matrixFile | cut -d "_" -f 1)
    val=$(python $countscript $matrixFile | grep "Result:" | cut -d " " -f 2)
    echo "$name,$val" >> $outFile
done


# Serio data
inDir="/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/beam"
outFile="/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/phylogenetic_information_to_cell_ratio_serio.csv"

echo "name,informative_muts_to_cell_ratio" > $outFile

matrixFiles=$(find $inDir -type f -name "temp_matrix.tsv")


for matrixFile in $matrixFiles; do
    name=$(echo $matrixFile | rev | cut -d "/" -f2-3 | rev | sed 's/\//_/g')
    val=$(python $countscript $matrixFile | grep "Result:" | cut -d " " -f 2)
    echo "$name,$val" >> $outFile
done

# Remove serio CP00 rows from the file after running the above script
sed -i '/CP00,/d' $outFile


# Simeonov data
inDir="/grid/siepel/home/staklins/stored_results/beam/latest_results/simeonov_2021_pancreatic_cancer_data/reformatted_raw_data/"
outFile="/grid/siepel/home/staklins/stored_results/beam/latest_results/simeonov_2021_pancreatic_cancer_data/phylogenetic_information_to_cell_ratio_simeonov.csv"

echo "name,informative_muts_to_cell_ratio" > $outFile

matrixFiles=$(find $inDir -type f -name "*matrix.csv")

for matrixFile in $matrixFiles; do
    name=$(basename $(dirname $matrixFile))_$(basename $matrixFile | cut -d "_" -f 1)
    val=$(python $countscript $matrixFile | grep "Result:" | cut -d " " -f 2)
    echo "$name,$val" >> $outFile
done

# Yang data
inDir="/grid/siepel/home/staklins/stored_data/crispr_barcode_related_data/yang_2022_real_data/processed_metadata/"
outFile="/grid/siepel/home/staklins/stored_results/beam/latest_results/yang_2022_lung_cancer_data/phylogenetic_information_to_cell_ratio_yang.csv"

echo "name,informative_muts_to_cell_ratio" > $outFile

matrixFiles=$(find $inDir -type f -name "*_successive_char_matrix_collapsed.tsv")

for matrixFile in $matrixFiles; do
    name=$(basename $matrixFile | cut -d "_" -f 1)-$(basename $matrixFile | cut -d "_" -f 2)-$(basename $matrixFile | cut -d "_" -f 3)
    val=$(python $countscript $matrixFile | grep "Result:" | cut -d " " -f 2)
    echo "$name,$val" >> $outFile
done
