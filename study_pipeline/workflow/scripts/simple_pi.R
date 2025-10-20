# Title     : pop gen stats
# Objective : Calculate Pi from gene alignment
# Created by: Evangelos A. Dimopoulos"
# Created on: 28/04/2025


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
aln_file <- args[1]
gene <- args[2]
out_file <- args[3]

# load libs
library(pegas)
library(ape)



# get an estimate of the nucleotide diversity (Pi) for each fasta file given
simple_pi <- function(aln_file, gene, out_file) {

    # read the fasta formatted alignment file
    dna <- read.dna(file = aln_file, format = "fasta", as.matrix = TRUE)

    # get pi for each alignment
    pi <- nuc.div(dna)

    # print to stdout the result for each fasta file
    print(paste(gene, round(pi, 3)), sep='    ')

    # save it to a file
    write(paste(gene, round(pi, 3), sep='    '),file=out_file,append=TRUE)

}

# run the function
simple_pi(aln_file, gene, out_file)