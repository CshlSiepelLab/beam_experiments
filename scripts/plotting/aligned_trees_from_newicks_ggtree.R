#!/usr/bin/env Rscript

library(ape)
library(ggtree)
library(gridExtra)
library(grid)
library(phytools)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# nwk1 <- args[1]
# nwk2 <- args[2]

nwk1 <- "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/laml/5k/28/laml_trees_no_origin.nwk"
nwk2 <- "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/28/sampled_tree_1_with_tissue_appended_to_name.nwk"

outfile <- sub("\\.nwk$", ".pdf", nwk2)

# Read in trees from provided Newick file paths
left_tree <- read.tree(nwk1)
right_tree <- read.tree(nwk2)

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
    tissues <- sapply(strsplit(tip_labels, "\\."), `[`, 1)
    tissues <- gsub("[0-9]", "", tissues)    # Collapse multiple samples from the same tissue into one tissue label
    tissues <- ifelse(grepl("R", tissues), "RL", tissues)   # Use coarse grained tissue labels for right lung
    return(tissues)
}

left_tree_tips <- data.frame(label = left_tree$tip.label)
left_tree_tips$tissue <- extract_tissue(left_tree_tips$label)

right_tree_tips <- data.frame(label = right_tree$tip.label)
right_tree_tips$tissue <- extract_tissue(right_tree_tips$label)

# Define a consistent color palette for tissues
tissues_all <- unique(c(left_tree_tips$tissue, right_tree_tips$tissue))

palette <- c(
    "LL" = "black",
    "M" = "red",
    "RL" = "blue",
    "Liv" = "green"
)
# For any tissues not listed above, assign a default color (e.g., "grey")
missing_tissues <- setdiff(tissues_all, names(palette))
if (length(missing_tissues) > 0) {
    palette <- c(palette, setNames(rep("grey", length(missing_tissues)), missing_tissues))
}

# To plot without the labels cutoff
max_x1 <- max(nodeHeights(left_tree)) + 10
max_x2 <- max(nodeHeights(right_tree)) + 10

# Plot left tree mirrored
p1 <- ggtree(left_tree) %<+% left_tree_tips +
  geom_tiplab(align = TRUE, linetype = 'dotted', hjust = 0, size = 0.75) +
  geom_tippoint(aes(color = tissue), size = 0.5) +
  scale_color_manual(values = palette) +
  theme_tree2() +
  theme(plot.margin = unit(c(1, 6, 1, 1), "lines"), legend.position = "none") +
  xlim(0, max_x1) +
  coord_cartesian(clip = 'off')

p2 <- ggtree(right_tree) %<+% right_tree_tips +
  geom_tiplab(align = TRUE, linetype = 'dotted', hjust = 1, size = 0.75) +
  geom_tippoint(aes(color = tissue), size = 0.5) +
  scale_color_manual(values = palette) +
  scale_x_reverse() +
  theme_tree2() +
  theme(plot.margin = unit(c(1, 1, 1, 6), "lines"), legend.position = "none") +
  xlim(max_x2, 0) +
  coord_cartesian(clip = 'off')

# Combine and save
pdf(outfile, width = 12, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()
