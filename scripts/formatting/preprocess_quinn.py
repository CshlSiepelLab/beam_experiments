import os
import sys
import pandas as pd
import cassiopeia as cas

from .matrix_utils import convert_matrix_to_row_successive_matrix


def convert_quinn_allele_table_to_successive_matrix(
    infile: str,
    lineage: int,
    outdir: str,
) -> None:
    """
    Converts a Quinn-style allele table to a successive character matrix and writes output files.
    This function reads an allele table from a TSV file, filters it for a specified lineage group,
    computes empirical indel priors, and converts the filtered allele table into both an original
    and a successive character matrix. It then writes the resulting matrices, mutation priors, and
    mutation dictionaries to the specified output directory.
    
    The pre-processing generally follows the guidelines outlined in the Cassiopeia documentation
    at https://cassiopeia-lineage.readthedocs.io/en/latest/notebooks/reconstruct.html to get the
    idel priors by sharing information across intBCs and then make a character matrix from the allele
    table by filtering low information sites.

    Args:
        infile (str): Path to the input allele table file (TSV format).
        lineage (int): The lineage group identifier to filter the allele table.
        outdir (str): Directory where output files will be written.

    Outputs:
        - <outdir>/<lineage>_mutation_priors.txt: Successive edit rates for each mutation code.
        - <outdir>/<lineage>_original_character_matrix.tsv: Original character matrix for the lineage group.
        - <outdir>/<lineage>_original_chracter_int_to_mutation_string_dict.txt: Mapping of character integers to mutation strings for the original matrix.
        - <outdir>/<lineage>_successive_character_matrix.tsv: Successive character matrix for the lineage group.
        - <outdir>/<lineage>_successive_int_to_mutation_string_dict.txt: Mapping of character integers to mutation strings for the successive matrix.

    """
    # read in the provided allele table
    allele_table = pd.read_csv(
        infile,
        sep="\t",
        usecols=[
            "cellBC",
            "intBC",
            "r1",
            "r2",
            "r3",
            "allele",
            "LineageGroup",
            "sampleID",
            "readCount",
            "UMI",
        ],
    )

    group = allele_table[allele_table["LineageGroup"] == lineage]

    # Get indel priors as per Cassiopeia docs - the cassiopeia
    # docs use the full allele table but group by lineage for 
    # counting whereas here we subset to the lineage first as
    # input, so these should be close
    indel_priors = cas.pp.compute_empirical_indel_priors(group)

    char_matrix_df, priors, mut_dict = cas.pp.convert_alleletable_to_character_matrix(
        group,
        missing_data_state="-1",
        allele_rep_thresh=0.95, # From Cassiopeia docs, but gives same number of sites as in the Quinn et al. deposited character matrices
        mutation_priors=indel_priors,
    )

    successive_matrix, new_mut_dict, successive_edit_rates = (
        convert_matrix_to_row_successive_matrix(char_matrix_df, mut_dict, indel_priors)
    )

    os.makedirs(outdir, exist_ok=True)

    # write successive edit rates to a file
    successive_edit_rates = dict(sorted(successive_edit_rates.items()))
    with open(f"{outdir}/{lineage}_mutation_priors.txt", "w") as f:
        f.write(f"mutation_code,rate\n")
        for key, value in successive_edit_rates.items():
            f.write(f"{key},{value}\n")

    # write each lineage group's successive matrix to its own file
    char_matrix_df.to_csv(
        f"{outdir}/{lineage}_original_character_matrix.tsv",
        sep="\t",
        index=True,
        header=True,
    )

    with open(
        f"{outdir}/{lineage}_original_chracter_int_to_mutation_string_dict.txt", "w"
    ) as f:
        f.write(f"site_num,char_int,mut_str\n")
        for key, value in mut_dict.items():
            for k, v in value.items():
                f.write(f"{key},{k},{v}\n")

    successive_matrix.to_csv(
        f"{outdir}/{lineage}_successive_character_matrix.tsv",
        sep="\t",
        index=True,
        header=True,
    )

    with open(
        f"{outdir}/{lineage}_successive_int_to_mutation_string_dict.txt", "w"
    ) as f:
        f.write(f"successive_char_int,mut_str\n")
        for key, value in new_mut_dict.items():
            f.write(f"{value},{key}\n")
            
def __main__():
    infile = sys.argv[1]
    lineage = int(sys.argv[2])
    outdir = sys.argv[3]
    
    convert_quinn_allele_table_to_successive_matrix(infile, lineage, outdir)