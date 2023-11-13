# Title     : root a tree
# Objective : Run host transition (switch) analysis
# Created by: Evangelos A. Dimopoulos"
# Created on: 23/06/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
outgroup <- args[2]
output <- args[3]
pruned <- args[4]

# load libs
library(ape)

# read the utility functions
source("scripts/utilities.R")

# treeWAS function
root_treefile <- function(tree_file, outgroup, output, pruned) {
    # read the tree and root it
    tree_r <- if (pruned=='pruned') read_tree(tree_file, outgroup) else root_tree(tree_file, outgroup)

    # write output table
    write.tree(tree_r, file=output)
}

# run the function
root_treefile(tree_file, outgroup, output, pruned)
