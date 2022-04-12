# Title     : heritability
# Objective : Calculate broad sense heritability (H^2)
# Created by: Evangelos A. Dimopoulos"
# Created on: 12/04/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
outgroup <- args[2]
aln_file_rec <- args[3]
aln_file_nrec <- args[4]
pop_metadata <- args[5]
out_table <- args[6]

# load libs
library(ade4)
library(adegenet)
library(ape)
library(poppr)

# read the utility functions
source("scripts/utilities.R")

# poppr function
heritability_calculate <- function(tree_file, outgroup, aln_file, pop_metadata) {

    # read tree
    tree_r_no_out <- read_tree(tree_file, outgroup)
    
    # read aln data from file:
    genind_obj <- read_aln(aln_file, tree_r_no_out)
    
    # read file with population metadata
    pop_meta <- population_host_metadata(pop_metadata)
    
    # load traits and define the population strata
    pop_trait <-
      pop_meta[pop_meta$sample %in% tree_r_no_out$tip.label, ]
    strata(genind_obj) <- pop_trait[, c("Trait", "Host", "Source")]
    setPop(genind_obj) <- ~ Trait
    
    # perform an amova analysis on the alignment
    amova_res <- poppr.amova(genind_obj, ~ Trait)
    
    # calculate n, m and H^2 based on http://ib.berkeley.edu/courses/ib162/Week4a.htm
    n = amova_res$results['Between samples', 'Df'] + 1
    m = (amova_res$results['Total', 'Df'] + 1) / n
    
    broad_h = ((amova_res$results['Between samples', 'Mean Sq'] - amova_res$results['Within samples', 'Mean Sq']) /
                 m) / amova_res$results['Total', 'Mean Sq']
    
    return(broad_h)
  }

# function to calculate heritability before and after masking recombination
heritability  <- function(tree_file, outgroup, aln_file_rec, aln_file_nrec, pop_metadata, out_table) {

    # get the H^2 for the alignments before we remove the recombination
    rec_h_broad <- heritability_calculate(tree_file, outgroup, aln_file_rec, pop_metadata)

    # get the H^2 for the alignments after we remove the recombination
    nrec_h_broad <- heritability_calculate(tree_file, outgroup, aln_file_nrec, pop_metadata)

    # combine the results
    rec_res <- c("Recombination", rec_h_broad)
    nrec_res <- c("No Recombination", nrec_h_broad)

    h_broad_res <- as.data.frame(rbind(rec_res, nrec_res))
    names(h_broad_res) <- c("Rec state", "Heritability (broad sense)")

    # write output table
    write.table(h_broad_res, file=out_table, row.names=FALSE, quote=FALSE, sep='\t')
}

# run the install function
heritability(tree_file, outgroup, aln_file_rec, aln_file_nrec, pop_metadata, out_table)