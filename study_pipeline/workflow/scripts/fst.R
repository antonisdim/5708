# Title     : treeWAS
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
out_table <- args[6]

# load libs
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

    return(c(nei_fstat['Fst'], wc_fstat$FST))
}

# function to calculate Fsts before and after masking recombination
fst  <- function(tree_file, outgroup, aln_file_rec, aln_file_nrec, pop_metadata, out_table) {

    # get the Fsts for the alignments before we remove the recombination
    rec_fst <- hierfstat_calculate(tree_file, outgroup, aln_file_rec, pop_metadata)

    # get the Fsts for the alignments after we remove the recombination
    nrec_fst <- hierfstat_calculate(tree_file, outgroup, aln_file_nrec, pop_metadata)

    # combine the results
    rec_res <- c("Recombination", rec_fst)
    nrec_res <- c("No Recombination", nrec_fst)

    fst_res <- as.data.frame(rbind(rec_res, nrec_res))
    names(fst_res) <- c("Rec state", "Nei's Fst", "WC's Fst")

    # write output table
    write.table(fst_res, file=out_table, row.names=FALSE, quote=FALSE, sep='\t')
}

# run the install function
fst(tree_file, outgroup, aln_file_rec, aln_file_nrec, pop_metadata, out_table)