#!/usr/bin/env python3

import scanpy as sc
import scipy.io
import pandas as pd
from Bio import Phylo
import os
from cassiopeia.data import CassiopeiaTree


input_matrix = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/data/quinn_2021_real_data/GSE161363/GSM4905335_matrix.5k.mtx"
input_barcodes = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/data/quinn_2021_real_data/GSE161363/GSM4905335_barcodes.5k.tsv"
input_genes = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/data/quinn_2021_real_data/GSE161363/GSM4905335_genes.5k.tsv"


# Load the matrix
matrix = scipy.io.mmread(input_matrix).T.tocsr()  # transpose if cells x genes is expected

# Load barcodes and gene names
barcodes = pd.read_csv(input_barcodes, header=None)[0].values
genes = pd.read_csv(input_genes, sep="\t", header=None)[0].values  # second column has gene names

# Create AnnData object
adata = sc.AnnData(X=matrix)
adata.obs_names = barcodes
adata.var_names = genes

# Label groups
group = []
for cell in adata.obs.index:
    if cell[0] != "L":
        group.append(cell[0])
    else:
        group.append(cell.split(".")[0])

adata.obs["group"] = group

# Preprocessing
#sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)
#sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)


# DEA
sc.tl.rank_genes_groups(adata, groupby="group", method='wilcoxon')  # or 't-test', 'logreg'
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

de_results = adata.uns["rank_genes_groups"]
significant_gene_sets = {}
for label in set(group):
    df_de = pd.DataFrame({
        "gene": de_results["names"][label],
        "log2FC": de_results["logfoldchanges"][label],
        "p_value": de_results["pvals"][label],
        "q_value": de_results["pvals_adj"][label],
        "signif": de_results["pvals_adj"][label] < 0.05
    })
    df_de.to_csv(f"./DEA_results_{label}.tsv", sep="\t", index=False)
    #print(sum(df_de["signif"].values))

    # Filter for significant genes (FDR < 0.05)
    significant_genes = set(
        df_de[(df_de["q_value"] < 0.05) & (df_de["log2FC"] > 0.5)]["gene"]
    )
    significant_gene_sets[label] = significant_genes
    print(f"{label}: {len(significant_genes)} significant genes.")

# Intersection across all groups
common_significant_genes = set.intersection(*significant_gene_sets.values())
print(f"{len(common_significant_genes)} significant genes in common.")
all_significant_genes = set.union(*significant_gene_sets.values())
print(f"{len(all_significant_genes)} significant genes in total.")

# Write to file
pd.Series(sorted(all_significant_genes)).to_csv("DEA_all_signif_genes.tsv", sep="\t", index=False, header=False)


# Count matrix for each tree
for i in range(24, 100):
    filename = str(i) + ".nwk"
    if filename in os.listdir("./LAML_tree/"):
        tree = CassiopeiaTree()
        tree.populate_tree("./LAML_tree/" + filename)

        # Select cells and genes
        cells_selected = [cell for cell in tree.leaves if cell in adata.obs_names]
        genes_selected = [gene for gene in all_significant_genes if gene in adata.var_names]
        
        # Subset to selected cells and genes
        subset = adata[cells_selected, genes_selected]
        sc.pp.filter_genes(subset, min_cells=len(tree.leaves)*0.1)
        
        # Convert sparse matrix to dense (if needed)
        count_matrix = pd.DataFrame(
            subset.X.toarray() if hasattr(subset.X, "toarray") else subset.X, 
            index=subset.obs_names, 
            columns=subset.var_names
        )
        count_matrix.to_csv(f"tree{i}_readcounts.tsv", sep="\t")
        print(i, count_matrix.shape)