library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(ape)
library(ggimage)

pdf(NULL)

# default colors taken from metient method for consistency in visualizations
DEFAULT_COLORS <- rep(c("#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0", "brown", "black", "darkgreen", "purple", "blue"), 3)



# treefile <- commandArgs(trailingOnly = TRUE)[1]
# primary_tissue <- commandArgs(trailingOnly = TRUE)[2]

treefile <- "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/yang_2022_real_data_8_22_24/metastabayes/3457_Apc_T4/combined.mcc.tree"
primary_tissue <- "T"

output_file <- sub("\\.tree$", "_ggtree.pdf", treefile)

beast_tree <- read.beast(treefile)

beast_tree_df <- as_tibble(beast_tree)
beast_tree_df <- mutate(beast_tree_df, label = ifelse(is.na(label), as.character(node), label))
beast_tree <- as.phylo(beast_tree_df)

locationset <- beast_tree_df$location.set
locationsetprob <- beast_tree_df$location.set.prob

location_order <- c(primary_tissue, sort(unique(beast_tree_df$location[beast_tree_df$location != primary_tissue])))
num_locations <- length(location_order)
colors <- setNames(c("black", DEFAULT_COLORS[1:length(location_order)]), location_order)

location_df <- beast_tree_df[,c("location.set", "location.set.prob")]

pie_df <- data.frame(matrix(ncol = length(location_order), nrow = 0))
colnames(pie_df) <- location_order

for (i in seq_len(nrow(location_df))) {
  locations <- unlist(location_df$location.set[i])
  probs <- as.numeric(unlist(location_df$location.set.prob[i]))
  pie_df[i, locations] <- probs
}

pie_df[is.na(pie_df)] <- 0

pie_df["node"] <- beast_tree_df["label"]

tip_labels <- beast_tree$tip.label
node_labels <- beast_tree$node.label

# Remove pie chart for fake nodes caused by converting polytomy to bifurcating tree with node that have branch length of 0
# Assumes the root is the first element of nodes since this is the efault behavior when reading in the tree
fake_nodes <- beast_tree_df$label[beast_tree_df$branch.length == 0]
fake_nodes <- fake_nodes[-1]
matching_rows <- pie_df$node %in% fake_nodes
pie_df[matching_rows, -which(names(pie_df) == "node")] <- 0

pies <- nodepie(pie_df, cols = 1:num_locations)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = colors))

node_pies <- pies[node_labels]
tip_pies <- pies[tip_labels]

# # Make tree ultrametric for visualization consistency
# beast_tree <- force.ultrametric(beast_tree)

max_width <- as.numeric(max(beast_tree_df$height))
x <- max_width * 1.5

p <- ggtree(beast_tree, aes(color=beast_tree_df$location), layout="rectangular") +
    geom_tippoint(aes(color = beast_tree_df$location), size=2) + 
    scale_color_manual(values = colors, breaks = location_order) +
    labs(color = "Tissue") +
    theme(legend.position = "right") +
    geom_inset(node_pies, width = 0.03)

ggsave(output_file, p)