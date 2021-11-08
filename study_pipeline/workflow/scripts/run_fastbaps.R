# Title     : run_fastbaps
# Objective : Cluster genomes with fastbaps
# Created by: Evangelos A. Dimopoulos"
# Created on: 08/11/2021


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
msa <- args[1]
cluster_file <- args[2]
threads <- as.integer(args[3])

#load libs
library(fastbaps)
library(ape)


run_fastbaps <- function(aln_path, out_path, threads) {
  # load snp msa
  sparse.data <- import_fasta_sparse_nt(aln_path)

  # apply baps prior
  sparse.data <- optimise_prior(sparse.data, type = "baps", n.cores = threads)

  # run fastbaps
  baps.hc <- fast_baps(sparse.data, n.cores = threads)

  # get clusters
  clusters <- best_baps_partition(sparse.data, as.phylo(baps.hc))

  # write output file
  write.table(clusters, out_path, col.names = FALSE, sep = ",", quote = FALSE)

}

# run the clustering function
run_fastbaps(msa, cluster_file, threads)
