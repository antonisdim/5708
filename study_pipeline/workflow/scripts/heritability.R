# Title     : heritability
# Objective : Calculate broad sense heritability (H^2)
# Created by: Evangelos A. Dimopoulos"
# Created on: 12/04/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
outgroup <- args[1]
aln_file <- args[2]
pop_metadata <- args[3]
out_table <- args[4]

# load libs
library(ade4)
library(adegenet)
library(ape)
library(poppr)

# read the utility functions
source("scripts/utilities.R")

# poppr function
heritability_calculate <- function(outgroup, aln_file, pop_metadata) {
    
    # read aln data from file:
    genind_obj <- read_aln(aln_file, outgroup, tree = FALSE)
    
    # read file with population metadata
    pop_meta <- population_host_metadata(pop_metadata)
    
    # load traits and define the population strata IN THE CORRECT ORDER - same order as in the genind object
    pop_trait <- pop_meta[pop_meta$sample %in% indNames(genind_obj[indNames(genind_obj) != outgroup]),]
    pop_trait <- pop_trait[match(indNames(genind_obj), pop_trait$sample),]
    strata(genind_obj) <- pop_trait[, c("Trait", "Host", "Source")]
    setPop(genind_obj) <- ~ Trait
    
    # perform an amova analysis on the alignment
    amova_res <- list()
    tryCatch({amova_res <- poppr.amova(genind_obj, ~ Trait, cutoff = 0.1)},
            error= function(e) {cat("ERROR :",conditionMessage(e), "\n",
            "H^2 cannot be calculated for that region. There are too many missing sites (>10%).",
            "\n")})

    if (length(amova_res) == 0) {
        broad_h <- 0
    } else {
        # calculate n, m and H^2 based on http://ib.berkeley.edu/courses/ib162/Week4a.htm
        n <- amova_res$results['Between samples', 'Df'] + 1
        m <- (amova_res$results['Total', 'Df'] + 1) / n
    
        broad_h <- ((amova_res$results['Between samples', 'Mean Sq'] - amova_res$results['Within samples', 'Mean Sq']) /
                    m) / amova_res$results['Total', 'Mean Sq']
    }

    return(broad_h)
  }

# function to calculate heritability before and after masking recombination
heritability  <- function(outgroup, aln_file, pop_metadata, out_table) {

    # get the H^2 for the alignments
    h_broad <- heritability_calculate(outgroup, aln_file, pop_metadata)

    # check if the alignment has recombinant sites masked
    rec_state <- if (grepl('nrec', aln_file)) "No Recombination" else "Recombination"

    # format the results
    res_vector <- c(rec_state, h_broad)
    h_broad_res <- as.data.frame(t(res_vector))
    names(h_broad_res) <- c("Rec state", "Heritability (broad sense)")

    # write output table
    write.table(h_broad_res, file=out_table, row.names=FALSE, quote=FALSE, sep='\t')
}

# run the install function
heritability(tree_file, outgroup, aln_file, pop_metadata, out_table)