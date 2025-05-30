#!/usr/bin/env Rscript

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(enrichplot)

# Read in the differential expression results
de_results <- read.delim("/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/gene_expression_analysis_rl_no_rl/rl_vs_norl_large_fc_all_genes.tsv", sep="\t")

# Filter for significant genes
threshold <- 0.05
sig_genes <- de_results[de_results$qval < threshold, ]

# Split genes into up and down regulated
up_genes <- sig_genes[sig_genes$log2fc > 0, "gene"]
down_genes <- sig_genes[sig_genes$log2fc < 0, "gene"]

# Create output directory if it doesn't exist
output_dir <- "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/gene_expression_analysis_rl_no_rl/"

# Function to run GO enrichment analysis
run_go_analysis <- function(gene_list, name) {
  # Run GO enrichment
  ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  ego_simple <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  return(ego_simple)
}

# Run analysis for up and down regulated genes
up_ego <- run_go_analysis(up_genes, "up_regulated")
down_ego <- run_go_analysis(down_genes, "down_regulated")

# Create various plots for up-regulated genes
up_dotplot <- dotplot(up_ego, 
                      showCategory = 10,
                      label_format = 70) + 
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

up_barplot <- barplot(up_ego, 
                     showCategory = 10,
                     label_format = 70) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

up_emapplot <- emapplot(pairwise_termsim(up_ego), 
                       showCategory = 10) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

up_cnetplot <- cnetplot(up_ego, 
                       showCategory = 5,
                       node_label = "gene") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# Create various plots for down-regulated genes
down_dotplot <- dotplot(down_ego, 
                       showCategory = 10,
                       label_format = 70) + 
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

down_barplot <- barplot(down_ego, 
                       showCategory = 10,
                       label_format = 70) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

down_emapplot <- emapplot(pairwise_termsim(down_ego), 
                         showCategory = 10) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

down_cnetplot <- cnetplot(down_ego, 
                         showCategory = 5,
                         node_label = "gene") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# Save all plots
# Up-regulated plots
ggsave(file.path(output_dir, "RL_enriched_GO_terms_dotplot.pdf"), 
       up_dotplot,
       width = 12, 
       height = 10, 
       dpi = 300,
       device = cairo_pdf)

ggsave(file.path(output_dir, "RL_enriched_GO_terms_barplot.pdf"), 
       up_barplot,
       width = 12, 
       height = 10, 
       dpi = 300,
       device = cairo_pdf)

ggsave(file.path(output_dir, "RL_enriched_GO_terms_emapplot.pdf"), 
       up_emapplot,
       width = 12, 
       height = 10, 
       dpi = 300,
       device = cairo_pdf)

ggsave(file.path(output_dir, "RL_enriched_GO_terms_cnetplot.pdf"), 
       up_cnetplot,
       width = 12, 
       height = 10, 
       dpi = 300,
       device = cairo_pdf)

# Down-regulated plots
ggsave(file.path(output_dir, "noRL_enriched_GO_terms_dotplot.pdf"), 
       down_dotplot,
       width = 12, 
       height = 10, 
       dpi = 300,
       device = cairo_pdf)

ggsave(file.path(output_dir, "noRL_enriched_GO_terms_barplot.pdf"), 
       down_barplot,
       width = 12, 
       height = 10, 
       dpi = 300,
       device = cairo_pdf)

ggsave(file.path(output_dir, "noRL_enriched_GO_terms_emapplot.pdf"), 
       down_emapplot,
       width = 12, 
       height = 10, 
       dpi = 300,
       device = cairo_pdf)

ggsave(file.path(output_dir, "noRL_enriched_GO_terms_cnetplot.pdf"), 
       down_cnetplot,
       width = 12, 
       height = 10, 
       dpi = 300,
       device = cairo_pdf)

# Save results to files
write.table(as.data.frame(up_ego), 
            file.path(output_dir, "RL_enriched_GO_terms.tsv"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(as.data.frame(down_ego), 
            file.path(output_dir, "noRL_enriched_GO_terms.tsv"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

