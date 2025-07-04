#!/usr/bin/env Rscript

library(ape)
library(ggtree)
library(gridExtra)
library(grid)
library(phytools)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggtreeExtra)
library(readr)
library(tidyr)
library(tibble)
library(ggnewscale)

# args <- commandArgs(trailingOnly = TRUE)
# nwk1 <- args[1] # LAML tree with branch lengths but no origin
# nwk2 <- args[2] # BEAM tree
# mutations <- args[3]

nwk1 <- "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/laml/5k/52/laml_trees_no_origin.nwk"
nwk2 <- "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/52/sampled_tree_3_with_tissue_appended_to_name.nwk"
mutations <- "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/successive_raw_data/5k/52_successive_character_matrix.tsv"

outfile <- sub("\\.nwk$", "_with_mutations.pdf", nwk2)

# Read in trees from provided Newick file paths
left_tree <- read.tree(nwk1)
right_tree <- read.tree(nwk2)

# Read in mutation data
mut_mat <- read_tsv(mutations)
mut_mat <- as.data.frame(mut_mat)
rownames(mut_mat) <- mut_mat[[1]]
mut_mat[[1]] <- NULL

# If no branch lengths, assign fake ones to allow layout
if (is.null(left_tree$edge.length)) {
  left_tree$edge.length <- rep(1, nrow(left_tree$edge))
}
if (is.null(right_tree$edge.length)) {
  right_tree$edge.length <- rep(1, nrow(right_tree$edge))
}

# Reorder tips in right_tree to match left_tree tip labels order
right_tree_ordered <- rotateConstr(right_tree, left_tree$tip.label)
right_tree <- right_tree_ordered

# Helper function to extract tissue from tip labels
extract_tissue <- function(tip_labels) {
    tissues <- ifelse(grepl("\\.", tip_labels),
                      sapply(strsplit(tip_labels, "\\."), `[`, 1),
                      "None")
    tissues <- gsub("'", "", tissues)  # Remove any single quote characters
    tissues <- gsub("[0-9]", "", tissues)    # Collapse multiple samples from the same tissue into one tissue label
    tissues <- ifelse(grepl("R", tissues), "RL", tissues)   # Use coarse grained tissue labels for right lung
    return(tissues)
}

left_tree_labels <- rbind(data.frame(label = left_tree$node.label), 
                            data.frame(label = left_tree$tip.label)
                            )
left_tree_labels$tissue <- extract_tissue(left_tree_labels$label)

right_tree_labels <- rbind(data.frame(label = right_tree$node.label), 
                            data.frame(label = right_tree$tip.label)
                            )
right_tree_labels$tissue <- extract_tissue(right_tree_labels$label)

# Define a consistent color palette for tissues
tissues_all <- unique(c(left_tree_labels$tissue, right_tree_labels$tissue))

# For quinn data tissue labels
palette <- c(
    "LL" = "black",
    "M" = "red",
    "RL" = "blue",
    "Liv" = "green",
    "None" = "grey"
)

# Align mutations to the right tree
mut_mat_left <- mut_mat[left_tree$tip.label, ] %>%
    rownames_to_column("label") %>%
    pivot_longer(-label, names_to = "mutation", values_to = "state")

# Align mutations to the right tree
mut_mat_right <- mut_mat[right_tree$tip.label, ] %>%
    rownames_to_column("label") %>%
    pivot_longer(-label, names_to = "mutation", values_to = "state")


mut_palette <- c("-1" = "grey", "0" = "white", "1" = "black")
unique_states <- unique(as.character(mut_mat_right$state))
new_states <- setdiff(unique_states, names(mut_palette))
if (length(new_states) > 0) {
  extra_colors <- grDevices::rainbow(length(new_states))
  names(extra_colors) <- new_states
  mut_palette <- c(mut_palette, extra_colors)
}

# Clean names of tips
clean_tip_name <- function(x) {
  x <- sub("^[^.]*\\.", "", x)   # Remove everything but the cell barcode
  x <- gsub("-1", "", x)         # Remove '-1'
  return(x)
}

left_tree$tip.label <- clean_tip_name(left_tree$tip.label)
right_tree$tip.label <- clean_tip_name(right_tree$tip.label)
mut_mat_left$label <- clean_tip_name(mut_mat_left$label)
mut_mat_right$label <- clean_tip_name(mut_mat_right$label)
left_tree_labels$label <- clean_tip_name(left_tree_labels$label)
right_tree_labels$label <- clean_tip_name(right_tree_labels$label)

# To plot without the labels cutoff
max_x1 <- max(nodeHeights(left_tree)) + 50
max_x2 <- max(nodeHeights(right_tree)) + 50

# Plot the trees
p1 <- ggtree(left_tree) %<+% left_tree_labels +
    geom_tree(aes(color = tissue)) +
    geom_tiplab(align = TRUE, linetype = 'dotted', linesize = 0, hjust = 0, size = 1, offset = 2) +
    geom_tippoint(aes(color = tissue), size = 0.5) +
    scale_color_manual(values = palette) +
    theme_tree2() +
    theme(plot.margin = unit(c(1, 2, 1, 1), "lines"), legend.position = "none") +
    xlim(0, max_x1) +
    coord_cartesian(ylim = c(-1, NA), clip = 'off')

p2 <- ggtree(right_tree) %<+% right_tree_labels +
    geom_tree(aes(color = tissue)) +
    geom_tiplab(align = TRUE, linetype = 'dotted', linesize = 0, hjust = 0, size = 1, offset = 2) +
    geom_tippoint(aes(color = tissue), size = 0.5) +
    scale_color_manual(values = palette) +
    theme_tree2() +
    theme(plot.margin = unit(c(1, 1, 1, 2), "lines"), legend.position = "none") +
    xlim(0, max_x2) +
    coord_cartesian(ylim = c(-1, NA), clip = 'off')

# Add the mutations
p1 <- p1 +
    geom_fruit(
      data = mut_mat_left,
      geom = geom_tile,
      mapping = aes(y = label, x = mutation, fill = as.factor(state)),
      offset = 0.2,
      pwidth = 0.5,
      color = "darkgrey",
      size = 0.01
    ) +
    scale_fill_manual(values = mut_palette, na.value = "white", drop = FALSE)

p2 <- p2 +
    geom_fruit(
      data = mut_mat_right,
      geom = geom_tile,
      mapping = aes(y = label, x = mutation, fill = as.factor(state)),
      offset = 0.2,
      pwidth = 0.5,
      color = "darkgrey",
      size = 0.01
    ) +
    scale_fill_manual(values = mut_palette, na.value = "white", drop = FALSE)

# Combine and save
pdf(outfile, width = 8, height = 8)
grid.arrange(p1, p2, ncol = 2)
dev.off()
