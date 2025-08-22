
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns


input_matrix = "/grid/siepel/home/staklins/stored_data/crispr_barcode_related_data/quinn_2021_real_data/GSE161363/GSM4905335_matrix.5k.mtx"
input_barcodes = "/grid/siepel/home/staklins/stored_data/crispr_barcode_related_data/quinn_2021_real_data/GSE161363/GSM4905335_barcodes.5k.tsv"
input_genes = "/grid/siepel/home/staklins/stored_data/crispr_barcode_related_data/quinn_2021_real_data/GSE161363/GSM4905335_genes.5k.tsv"
input_meta = "/grid/siepel/home/staklins/stored_data/crispr_barcode_related_data/quinn_2021_real_data/GSE161363/GSM4905335_meta.5k.tsv"

outfile = "/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data/rl_vs_norl_sig_gene_log2fc_barplot.pdf"

# Load count matrix
adata = sc.read(input_matrix, cache=True).T	# cells as rows, genes as columns

# Load in gene names and add to adata
genes = pd.read_csv(input_genes, header=None, sep='\t')
adata.var_names = genes[1]
adata.var_names = [x.upper() for x in adata.var_names]
adata.var_names_make_unique()

# Load in cell names and add to adata
barcodes = pd.read_csv(input_barcodes, header=None)
adata.obs_names = barcodes[0]

# Load in the cell metadata, keep only cells in adata that have metadata, and add metadata to adata
meta = pd.read_csv(input_meta, sep='\t', index_col=0)
ts_rna_overlap = np.intersect1d(adata.obs_names, meta.index)
meta = meta.loc[ts_rna_overlap]
adata = adata[ts_rna_overlap,:]
adata.obs = adata.obs.merge(meta, left_index = True, right_index=True, how="left")

# Filter cells with > 20% mitochondrial gene expression
mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
mito_frac = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata = adata[mito_frac <= 0.20, :]

# Filter cells with few genes expressed
sc.pp.filter_cells(adata, min_genes=100)

# Filter out genes that are not expressed in at least 1% of cells
sc.pp.filter_genes(adata, min_cells=0.01*adata.shape[0])

# Scale the library size to counts per million based on median total counts per cell
sc.pp.normalize_total(adata)

# Log transform the data
sc.pp.log1p(adata)

# Subset to only the LL data
adata_ll = adata[adata.obs.apply(lambda x: x['sampleID'] == "LL", axis=1),:]

# Set the CPs/groups that chose the RL model and no RL model based on external hypothesis tests
groups_RL_model = np.array(['32', '34', '43', '37', '26', '47', '40', '30', '79', '63', '36', '58', '70', '56', '57', '55', '61', '66', '84', '62', '80', '83', '74', '67', '86', '98', '52', '90', '92', '91', '95', '68'], dtype=int)
groups_no_RL_model = np.array(['42', '24', '35', '28', '45', '44', '51', '82', '64', '59', '72', '60', '89', '71', '73', '99', '97', '77', '76', '100', '54'], dtype=int)

# Subset to only the groups in either model and label cells with groups
adata_ll = adata_ll[adata_ll.obs.apply(lambda x: x['LineageGroup'] in groups_RL_model or x['LineageGroup'] in groups_no_RL_model, axis=1),:]
adata_ll.obs["model"] = adata_ll.obs.apply(lambda x: "RL" if x['LineageGroup'] in groups_RL_model else "noRL", axis=1)
adata_ll.obs["model"] = adata_ll.obs["model"].astype('category')

# Remove mitochondrial and ribosomal genes
remove_genes = [name for name in adata_ll.var_names if name.startswith('MT-') or name.startswith('RPL') or name.startswith('RPS')]
adata_ll = adata_ll[:, [gene for gene in adata_ll.var_names if gene not in remove_genes]].copy()

# Only keep genes in >20% of cells in at least one group
min_frac = 0.2
groupby_col = "model"
mask = np.any([
    (adata_ll[adata_ll.obs[groupby_col] == grp].X > 0).mean(axis=0) >= min_frac
    for grp in adata_ll.obs[groupby_col].unique()
], axis=0).flatten()
adata_ll = adata_ll[:, mask].copy()

# Run the DEA analysis
sc.tl.rank_genes_groups(adata = adata_ll, 
                        groupby = "model",
                        groups = ["RL"],
                        reference = "noRL",
                        method = "wilcoxon",
                        only_positive = False)

result = adata_ll.uns["rank_genes_groups"]

# Create a DataFrame for RL vs noRL results
rl_vs_norl = pd.DataFrame({
	'gene': result['names']['RL'],
	'logfc': result['logfoldchanges']['RL'],
	'log2fc': np.log2(np.exp(result['logfoldchanges']['RL'])),  # convert lnFC → log2FC
	'pval': result['pvals']['RL'],
	'qval': result['pvals_adj']['RL']
}).sort_values(
	by='qval',
	ascending=True
)

# Filter for significant genes (q-value < 0.05 and |log2FC| > log2(1.5))
rl_vs_norl_sig = rl_vs_norl[(rl_vs_norl['qval'] < 0.05) & (np.abs(rl_vs_norl['log2fc']) > np.log2(1.5))]
rl_vs_norl_sig = rl_vs_norl_sig.sort_values(by='log2fc', ascending=False)

# Write the significant genes to a file
rl_vs_norl_sig.to_csv(outfile.replace(".pdf", ".tsv"), sep="\t", index=False)

# Make a two sided barplot of log2fc of significant genes
fig = plt.figure(figsize=(6,10))
colors = ['#d62728' if x < 0 else '#1f77b4' for x in rl_vs_norl_sig['log2fc']]
sns.barplot(data=rl_vs_norl_sig, y='gene', x='log2fc', palette=colors)
plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
plt.ylabel('Genes')
plt.xlabel('Log2 Fold Change of RL vs noRL')
plt.title('Log2 Fold Change of Significant Genes')
plt.tight_layout()
plt.savefig(outfile, bbox_inches='tight')
plt.close()
