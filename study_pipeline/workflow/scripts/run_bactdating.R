# Title     : BactDating
# Objective : Date the nodes of the tree
# Created by: Evangelos A. Dimopoulos"
# Created on: 06/11/2023


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
treefile <- args[1]
outgroup <- args[2]
genome_length <- args[3]
pop_metadata <- args[4]
out_tree <- args[5]

# load libs
library(ape)
library(BactDating)
library(coda)
library(treeio)

# read the utility functions
source("scripts/utilities.R")

# define the function that runs bact dating
run_bactdating <- function(treefile, outgroup, genome_length, pop_metadata, out_tree) {

  # read the tree
  tree_obj <- read_tree(treefile, outgroup)

  # convert the branch lengths to mutations per genome (per year)
  tree_obj$edge.length <- tree_obj$edge.length * as.numeric(genome_length)

  # read file with population metadata
  pop_meta <- population_host_metadata(pop_metadata)

  # load traits and define the population strata IN THE CORRECT ORDER - same order as in the genind object
  pop_trait <- pop_meta[pop_meta$sample %in% tree_obj$tip.label,]
  pop_trait <- pop_trait[match(tree_obj$tip.label, pop_trait$sample),]

  # get the sampling dates
  collection_year <- as.numeric(pop_trait$Collection_Year)

  # run bactdating
  res <- bactdate(tree_obj, collection_year, nbIts=5e+7, showProgress=TRUE, updateRoot=FALSE)

  # convert bactdating results to a coda MCMC df
  mcmc <- as.mcmc.resBactDating(res, burnin = 0.1)
  ess <- effectiveSize(mcmc)

  print(res)
  print(ess)

  # convert the bactdating tree into a beast format
  tree_data <- as.treedata.resBactDating(res)
  beastlike_obj <- methods::new('treedata', phylo = tree_data[[1]],
                                data = dplyr::as_tibble(as.data.frame(tree_data[[2]])))
  write.beast(beastlike_obj, out_tree)
}

# run the install function
run_bactdating(treefile, outgroup, genome_length, pop_metadata, out_tree)