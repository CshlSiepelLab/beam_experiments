#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scanpy as sc
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns


def compute_log2fc(gr, bg, gene, adata):
	
	# g1_filt = adata.obs.apply(lambda x: x[groupby_var] == g, axis=1)
	# bg_filt = adata.obs.apply(lambda x: x[groupby_var] == bg, axis=1)
	
	gene_ii = np.where(adata.raw.var_names == gene)[0][0]
	
	exp_g = np.mean(gr[:,gene_ii]) + 0.01
	exp_bg = np.mean(bg[:,gene_ii]) + 0.01
	
	return np.log2(exp_g / exp_bg)


def create_DE_df(counts, groupby_var, gr, bg, result, method='ttest'):
	
	g_filt = counts.obs.apply(lambda x: x[groupby_var] ==  gr, axis=1).values
	bg_filt = counts.obs.apply(lambda x: x[groupby_var] == bg, axis=1).values
	
	print(counts.X.shape, len(g_filt), g_filt[:10])
	gdata = counts.X[g_filt, :]
	bgdata = counts.X[bg_filt, :]
	
	log2fc = {}
	adj_pvalues = {}
	scores = {}

	if method == 'logreg':
		for gene, score in zip(result['names'][gr], result['scores'][gr]):
			scores[gene] = score
			log2fc[gene] = compute_log2fc(gdata, bgdata, gene, counts)
		
		de_df = pd.DataFrame.from_dict(scores, orient='index', columns=['scores'])
		de_df['gene'] = de_df.index
		de_df['log2fc'] = de_df.index.map(log2fc)
		de_df.index = range(de_df.shape[0])
		return de_df
    	
	for gene, qval, fc in zip(result['names'][gr], result['pvals_adj'][gr], result['logfoldchanges'][gr]):
		#scores[gene] = score
		#log2fc[gene] = compute_log2fc(gdata, bgdata, gene)
		
		adj_pvalues[gene] = qval
		log2fc[gene] = compute_log2fc(gdata, bgdata, gene, counts)
		
	de_df = pd.DataFrame.from_dict(adj_pvalues, orient='index', columns=['qval'])
	de_df['gene'] = de_df.index
	de_df['log2fc'] = de_df.index.map(log2fc)
	de_df.index = range(de_df.shape[0])
	return de_df



input_matrix = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/data/quinn_2021_real_data/GSE161363/GSM4905335_matrix.5k.mtx"
input_barcodes = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/data/quinn_2021_real_data/GSE161363/GSM4905335_barcodes.5k.tsv"
input_genes = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/data/quinn_2021_real_data/GSE161363/GSM4905335_genes.5k.tsv"
input_meta = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/data/quinn_2021_real_data/GSE161363/GSM4905335_meta.5k.tsv"

# Load the count matrix and transpose it to get cells as rows and genes as columns
adata = sc.read(input_matrix, cache=True).T
# Load gene annotations from the genes file
genes = pd.read_csv(input_genes, header=None, sep='\t')
# Load cell barcodes from the barcodes file
barcodes = pd.read_csv(input_barcodes, header=None)

# Set gene names as variable names (columns) in the AnnData object
adata.var_names = genes[1]
# Store the original gene IDs as an annotation in the var dataframe
adata.var['gene_ids'] = genes[0]
# Set cell barcodes as observation names (rows) in the AnnData object
adata.obs_names = barcodes[0]
# Ensure all gene names are unique by appending numbers if necessary
adata.var_names_make_unique()

# Convert all gene names to uppercase for consistency
adata.var_names = [x.upper() for x in adata.var_names]

# Read in the metadata containing cell annotations and scores
meta = pd.read_csv(input_meta, sep='\t', index_col=0)

# Find cells that are present in both the expression data and metadata
ts_rna_overlap = np.intersect1d(adata.obs_names, meta.index)
# Filter metadata to only include cells present in expression data
meta = meta.loc[ts_rna_overlap]
# Filter expression data to only include cells present in metadata
adata = adata[ts_rna_overlap,:]

# Identify mitochondrial genes by finding genes that start with 'MT-'
mito_genes = [name for name in adata.var_names if name.startswith('MT-')]

# Calculate percentage of mitochondrial reads per cell
# This is a common quality control metric - high mitochondrial content often indicates cell stress or death
adata.obs['percent_mito'] = np.sum(adata[:,mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1


# Merge cell metadata with the AnnData object
# This adds all columns from meta to adata.obs, matching cells by their index
adata.obs = adata.obs.merge(meta, left_index = True, right_index=True, how="left")

# Calculate median library size across all cells to use as scaling factor
# Normalize cells to have the same total count (scale_factor)
# This helps account for differences in library size between cells
scale_factor = np.median(np.array(adata.X.sum(axis=1)))
sc.pp.normalize_per_cell(adata, counts_per_cell_after = scale_factor)

# Filter out cells with high mitochondrial content (>20%)
# High mitochondrial content often indicates cell stress or death
adata = adata[adata.obs.percent_mito <= 0.20, :]

# Filter out genes that are not expressed in at least 1% of cells
# This removes very rare genes that might be technical artifacts
thresh = 0.01*adata.shape[0]
sc.pp.filter_genes(adata, min_cells=thresh)

# Store the raw data before normalization for later use
adata.raw = adata

# Log transform the data (log1p = log(x + 1))
# This helps normalize the data distribution for downstream analysis
sc.pp.log1p(adata)

# Subset to only the LL data
adata_ll = adata[adata.obs.apply(lambda x: x['sampleID'] in ["LL"], axis=1),:]

# Set the CPs/groups that chose the RL model and no RL model based on external hypothesis tests

# # CPs from the threshold 3
# groups_RL_model = ['34', '43', '37', '47', '40', '30', '36', '70', '57', '60', '66', '62', '67', '74', '86', '71', '54', '24']
# groups_no_RL_model = ['42', '35', '28', '45', '44', '51', '79', '59', '64', '80', '82']

# CPs from the threshold 1.1
groups_RL_model = [
    "24", "26", "30", "32", "34", "36", "37", "40", "43", "47",
    "54", "57", "60", "62", "66", "67", "70", "71", "74", "84",
    "86", "91", "92", "98", "99"
]
groups_no_RL_model = [
    "28", "29", "35", "42", "44", "45", "46", "48", "49", "51", 
    "52", "55", "56", "58", "59", "61", "63", "64", "68", "72", 
    "73", "76", "77", "78", "79", "80", "82", "83", "85", "89", 
    "90", "94", "95", "96", "97", "100"
]

# make model lists int dtype
groups_RL_model = np.array(groups_RL_model, dtype=int)
groups_no_RL_model = np.array(groups_no_RL_model, dtype=int)

# Subset to only the groups in either model
adata_ll = adata_ll[adata_ll.obs.apply(lambda x: x['LineageGroup'] in groups_RL_model or x['LineageGroup'] in groups_no_RL_model, axis=1),:]

# Add a column to the adata_ll object to indicate which model each cell belongs to
adata_ll.obs["model"] = adata_ll.obs.apply(lambda x: "RL" if x['LineageGroup'] in groups_RL_model else "noRL", axis=1)

print("Number of cells in each group:")
print(adata_ll.obs["model"].value_counts())

# Convert model column to categorical type for differential expression analysis
adata_ll.obs["model"] = adata_ll.obs["model"].astype('category')

# Run the DEA analysis
sc.tl.rank_genes_groups(adata_ll, "model", groups=["RL"], reference="noRL", method="wilcoxon", use_raw = True, n_genes = len(adata_ll.var_names), only_positive = False)

# Get the results from the differential expression analysis
result = adata_ll.uns["rank_genes_groups"]
groups = result['names'].dtype.names

# Create differential expression dataframes for RL vs noRL comparison
rl_vs_norl = create_DE_df(adata_ll, "model", "RL", "noRL", result, method='wilcoxon')

# Filter for significant genes (q-value < 0.05 and |log2FC| > log2(1.5))
rl_vs_norl_sig = rl_vs_norl[(rl_vs_norl['qval'] < 0.05) & (np.abs(rl_vs_norl['log2fc']) > np.log2(1.5))]
rl_vs_norl_sig = rl_vs_norl_sig.sort_values(by='log2fc', ascending=False)

# Write the significant genes to a file
rl_vs_norl_sig.to_csv("rl_vs_norl_sig_genes.tsv", sep="\t", index=False)

# For the significant genes, get the proportion of cells that express the gene in each model and write to a file
prop_expr = []
for gene in rl_vs_norl_sig['gene']:
    # Get cells for each model
    rl_cells = adata_ll[adata_ll.obs['model'] == 'RL', gene]
    norl_cells = adata_ll[adata_ll.obs['model'] == 'noRL', gene]
    
    # Calculate percentage of cells with non-zero expression
    rl_prop = (rl_cells.X > 0).sum() / len(rl_cells)
    norl_prop = (norl_cells.X > 0).sum() / len(norl_cells)
    
    prop_expr.append({
        'gene': gene,
        'RL_prop': rl_prop,
        'noRL_prop': norl_prop,
    })

prop_expr_df = pd.DataFrame(prop_expr)
prop_expr_df.to_csv("rl_vs_norl_gene_expression_proportions.tsv", sep="\t", index=False)

# Get genes with large fold changes regardless of significance
rl_vs_norl_fc = rl_vs_norl[(np.abs(rl_vs_norl['log2fc']) > np.log2(1.5))]
rl_vs_norl_fc = rl_vs_norl_fc.sort_values(by='log2fc', ascending=False)

# Write the large fold change genes to a file
rl_vs_norl_fc.to_csv("rl_vs_norl_large_fc_all_genes.tsv", sep="\t", index=False)

# Get the significant genes
significant_genes = rl_vs_norl_sig['gene'].tolist()

# # Create matrix plot of significant genes
# fig = plt.figure(figsize=(5,6))
# ax = sc.pl.matrixplot(adata_ll, significant_genes, groupby='model', dendrogram=False,
#                 use_raw=False, log=True, vmin = -2, vmax=2,
#                 swap_axes=True, standard_scale='var', show=False)
# plt.savefig("rl_vs_norl_sig_gene_expression_matrix.pdf", bbox_inches='tight')
# plt.close()

# Make a two sided barplot of log2fc of significant genes
fig = plt.figure(figsize=(6,10))
# Create color palette based on log2fc values using muted colors
colors = ['#d62728' if x < 0 else '#1f77b4' for x in rl_vs_norl_sig['log2fc']]  # muted red and blue
sns.barplot(data=rl_vs_norl_sig, y='gene', x='log2fc', palette=colors)
plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
plt.ylabel('Genes')
plt.xlabel('Log2 Fold Change of RL vs noRL')
plt.title('Log2 Fold Change of Significant Genes')
plt.tight_layout()
plt.savefig("rl_vs_norl_sig_gene_log2fc_barplot.pdf", bbox_inches='tight')
plt.close()

# # Create dotplot of significant genes
# sc.set_figure_params(dpi=80, color_map='viridis')
# ax = sc.pl.dotplot(adata_ll, significant_genes, groupby='model',
#               figsize=(12,4),
#               standard_scale='var')


# # get the mean expression of a specific gene for each model
# gene = "KRT17"
# rl_mean = adata_ll[adata_ll.obs['model'] == 'RL', gene].X.mean(axis=0)
# norl_mean = adata_ll[adata_ll.obs['model'] == 'noRL', gene].X.mean(axis=0)

# print(rl_mean, norl_mean)


