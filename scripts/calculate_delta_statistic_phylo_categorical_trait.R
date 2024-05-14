########################################################################################
# Rui Borges, João Paulo Machado, Cidália Gomes, Ana Paula Rocha, Agostinho Antunes; Measuring phylogenetic signal between categorical traits and phylogenies, Bioinformatics, https://doi.org/10.1093/bioinformatics/bty800
# https://github.com/mrborges23/delta_statistic
########################################################################################

library(ape)
library(dplyr)
source("~/delta_statistic/code.R")

# # inputs
# newickfile <- commandArgs(trailingOnly = TRUE)[1]
# traitfile <- commandArgs(trailingOnly = TRUE)[2]

# for testing
# newickfile <- "results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/mS/1983/cell_tree_seed1983.nwk"
newickfile <- "results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/mS/11455/cell_tree_seed519673969.nwk"
traitfile <- "results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/mS/11455/cell_tree_seed519673969.labeling"

# Read the tree from a Newick file
tree <- read.tree(newickfile)

# check if tree is bifurcating
if (!is.binary(tree)) {
    print("The tree is not bifurcating. Converting to bifurcating topology.")
    tree <- multi2di(tree)
}

# replace 0 edge lengths with the minimum positive edge length
tree$edge.length[tree$edge.length <= 0] <- min(tree$edge.length[tree$edge.length > 0])

# Read trait data from a TSV file
trait_data <- read.table(traitfile, header=FALSE, sep=" ", col.names = c("tip_label", "trait"), colClasses = c("character", "character"))

# get tip names
tip_names <- tree$tip.label

# drop non tip names from trait df
trait_data <- trait_data[trait_data$tip_label %in% tip_names, ]

# add na for tips not in trait already
missing_tips <- setdiff(tip_names, trait_data$tip_label)
if (length(missing_tips) != 0) {
  print(paste("Missing tips: ", missing_tips))
  trait_data <- rbind(trait_data, data.frame(tip_label = missing_tips, trait = "undefined"))
}

# sort trait_data by tip_label to match tip_names order
trait_data <- trait_data %>%
  arrange(factor(tip_label, levels = tip_names))

# check if order is the same for trait_data$tip_label and tip_names
if (!identical(trait_data$tip_label, tip_names)) {
    stop("The order of trait_data$tip_label and tip_names is not the same.")
}

trait_data$trait <- as.numeric(as.factor(trait_data$trait))
trait <- trait_data$trait


# # used to downsample for testing with ape to check label matching
# num_downsample <- 10
# if(length(tree$tip.label) > num_downsample) {
#   # Get the names of the tips to keep
#   tips_to_keep <- sample(tree$tip.label, num_downsample)

#   # Get the names of the tips to drop
#   tips_to_drop <- setdiff(tree$tip.label, tips_to_keep)

#   # Drop the tips
#   tree <- drop.tip(tree, tips_to_drop)
# }
# trait <- trait_data[trait_data$tip_label %in% tree$tip.label, ]$trait


# required bug to remove tree$node.label for ace() from R package ape to run properly for discrete traits???
tree$node.label <- rep(NULL, length(tree$node.label))

# calculate the delta statistic
lambda0 <- 0.1   #rate parameter of the proposal 
se      <- 0.5   #standard deviation of the proposal
sim     <- 10000 #number of iterations
thin    <- 10    #we kept only each 10th iterate 
burn    <- 100   #100 iterates are burned-in
deltaA <- delta(trait,tree,lambda0,se,sim,thin,burn)
print(paste("Delta statistic: ", deltaA))

# # use below to plot the distribution of delta values under the null hypothesis
# num_reps <- 10
# random_delta <- rep(NA,num_reps)
# for (i in 1:num_reps){
#   rtrait <- sample(trait)
#   random_delta[i] <- delta(rtrait,tree,lambda0,se,sim,thin,burn)
# }
# p_value <- sum(random_delta>deltaA)/length(random_delta)
# boxplot(random_delta)
# abline(h=deltaA,col="red")