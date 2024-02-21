library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(ape)
library(ggimage)

treefile <- "beast_gundem_2015_2_21_24/A10_sym/tissue_tree_with_trait.tree"
primary_tissue <- "prostate" 

beast_tree <- read.beast(treefile)

beast_tree_df <- as_tibble(beast_tree)
beast_tree_df <- mutate(beast_tree_df, label = ifelse(is.na(label), as.character(node), label))
beast_tree <- as.phylo(beast_tree_df)

locationset <- beast_tree@data$location.set
locationsetprob <- beast_tree@data$location.set.prob

location_order <- c(primary_tissue, unique(beast_tree@data$location[beast_tree@data$location != primary_tissue]))
num_locations <- length(location_order)
colors <- setNames(palette()[1:length(location_order)], location_order)

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
fake_nodes <- beast_tree_df$label[beast_tree_df$branch.length == 0]
# Assumes the root is the first element of nodes since this is the efault behavior when reading in the tree
fake_nodes <- fake_nodes[-1]
matching_rows <- pie_df$node %in% fake_nodes
pie_df[matching_rows, -which(names(pie_df) == "node")] <- 0

pies <- nodepie(pie_df, cols = 1:num_locations)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = colors))

node_pies <- pies[node_labels]
tip_pies <- pies[tip_labels]

p <- ggtree(beast_tree) +
    geom_tiplab(hjust = -0.15) +
    geom_tippoint(aes(color = beast_tree_df$location), size=5) + 
    scale_color_manual(values = colors) +
    labs(color = "Tissue") +
    theme(legend.position = "right")

inset(p, node_pies, width=.08, height=.05)
