# Title     : fst
# Objective : Run Nei's (1982) and Weir and Cockerham's Fst
# Created by: Evangelos A. Dimopoulos"
# Created on: 08/04/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
outgroup <- args[2]
aln_file_rec <- args[3]
aln_file_nrec <- args[4]
pop_metadata <- args[5]
sum_out_table <- args[6]
pair_wc_fst_out_table <- args[7]
pair_nei_dist_out_table <- args[8]

# load libs
library(ade4)
library(adegenet)
library(ape)
library(hierfstat)

# read the utility functions
source("scripts/utilities.R")

# hierfstat function
hierfstat_calculate <- function(tree_file, outgroup, aln_file, pop_metadata) {

    # read tree
    tree_r_no_out <- read_tree(tree_file, outgroup)

    # read aln data from file:
    genind_obj <- read_aln(aln_file, tree_r_no_out)

    # read file with population metadata
    pop_meta <- population_host_metadata(pop_metadata)

    # load traits
    pop_trait <- pop_meta[pop_meta$sample %in% tree_r_no_out$tip.label,]
    phen_cat <- setNames(pop_trait$Trait, pop_trait$sample)

    # convert genind object to hierfstat input
    hierfstat_input <- genind2hierfstat(genind_obj, pop=phen_cat)

    # calculate Nei's Fst - based on genotype frequencies
    nei_fstat <- basic.stats(hierfstat_input, diploid=FALSE,digits=2)$overall

    # calculate WC's Fst - based on ANOVA and DNA sequence data
    wc_fstat <- wc(hierfstat_input, diploid=FALSE)

    # calculate WC's pairwise Fst - we prefer this to Nei's as the latter overestimates structure
    wc_pair_fstat <- pairwise.WCfst(hierfstat_input, diploid=FALSE)

    # calculate Nei's Ds - the 1984 formulation also known as Da
    nei_dist <- genet.dist(hierfstat_input, diploid=FALSE, method="Da")

    return(list(c(nei_fstat['Fst'], wc_fstat$FST), wc_pair_fstat, nei_dist))
}

# function to calculate Fsts before and after masking recombination
distances  <- function(tree_file, outgroup, aln_file_rec, aln_file_nrec, pop_metadata, sum_out_table,
    pair_wc_fst_out_table, pair_nei_dist_out_table) {

    # get the Fsts for the alignments before we remove the recombination
    rec_dist <- hierfstat_calculate(tree_file, outgroup, aln_file_rec, pop_metadata)

    # get the Fsts for the alignments after we remove the recombination
    nrec_dist <- hierfstat_calculate(tree_file, outgroup, aln_file_nrec, pop_metadata)

    # combine the results
    rec_res <- c("Recombination", rec_dist[[1]])
    nrec_res <- c("No Recombination", nrec_dist[[1]])

    fst_res <- as.data.frame(rbind(rec_res, nrec_res))
    names(fst_res) <- c("Rec state", "Nei's Fst", "WC's Fst")

    # write Fst summary output table
    write.table(fst_res, file=sum_out_table, row.names=FALSE, quote=FALSE, sep='\t')

    # write the pairwise WC Fst
    cat("WC pairwise Fst - with recombination", file=pair_wc_fst_out_table, sep="\n")
    write.table(as.matrix(rec_dist[[2]]), file=pair_wc_fst_out_table, quote=FALSE, row.names=TRUE, sep="\t",
        append=TRUE )
    cat("WC pairwise Fst - no recombination", file=pair_wc_fst_out_table, sep="\n", append=TRUE)
    write.table(as.matrix(nrec_dist[[2]]), file=pair_wc_fst_out_table, quote=FALSE, row.names=TRUE, sep="\t",
        append=TRUE )

    # write Nei's Ds output table
    cat("Nei's pairwise genetic distances - with recombination", file=pair_nei_dist_out_table, sep="\n")
    write.table(as.matrix(rec_dist[[3]]), file=pair_nei_dist_out_table, quote=FALSE, row.names=TRUE, sep="\t",
        append=TRUE )
    cat("Nei's pairwise genetic distances - no recombination", file=pair_nei_dist_out_table, sep="\n", append=TRUE)
    write.table(as.matrix(nrec_dist[[3]]), file=pair_nei_dist_out_table, quote=FALSE, row.names=TRUE, sep="\t",
        append=TRUE )
}

# run the install function
distances(tree_file, outgroup, aln_file_rec, aln_file_nrec, pop_metadata, sum_out_table,
    pair_wc_fst_out_table, pair_nei_dist_out_table)