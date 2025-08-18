
# Serio
indir="/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/beam"
outfile="/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data/cps_with_one_tissue_to_exclude.txt"

echo "mouse,cp" > $outfile

for file in $(find $indir -type f -name "expanded_clones_tissues.tsv"); do
    tissues=$(tail -n +2 $file | cut -f 2 | sort | uniq | wc -l)
    if [ "$tissues" -eq 1 ]; then
        mouse=$(basename $(dirname $(dirname $file)))
        cp=$(basename $(dirname $file))
        echo "$mouse,$cp" >> $outfile
    fi
done

# Quinn
indir="/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/beam"
outfile="/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/cps_with_one_tissue_to_exclude.txt"

echo "mouse,cp" > $outfile

for file in $(find $indir -type f -name "*tip_tissues.csv"); do
    tissues=$(cut -d',' -f 2 $file | sort | uniq | wc -l)
    if [ "$tissues" -eq 1 ]; then
        mouse=$(basename $(dirname $(dirname $file)))
        cp=$(basename $(dirname $file))
        echo "$mouse,$cp" >> $outfile
    fi
done