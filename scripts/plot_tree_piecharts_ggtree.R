library(ggtree)
library(treeio)
library(ggplot2)

treefile <- "gundem_a10/tissue_tree_with_trait.tree"
primary_tissue <- "prostate" 

beast_tree <- read.beast(treefile)

x <- c(primary_tissue, unique(beast_tree@data$location[beast_tree@data$location != primary_tissue]))

cols <- setNames(palette()[1:length(x)], x)

p <- ggtree(beast_tree) +
    geom_tiplab()

pies <- nodepie(beast_tree@data, cols = "location.set.prob")
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols))

p2 <- p + geom_inset(pies, width = .1, height = .1)
