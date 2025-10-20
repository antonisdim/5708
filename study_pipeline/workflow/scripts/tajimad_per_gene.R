# Title     : pop gen stats
# Objective : Run Nei's (1982) and Weir and Cockerham's Fst, Tajima's D, Nei's Dxy, and broad sense heritability H^2
# Created by: Evangelos A. Dimopoulos"
# Created on: 29/07/2029


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
aln_rec_chr <- args[1]
outgroup <- args[2]
pop_metadata <- args[3]
coords_file <- args[4]
output_table <- args[5]

# load libs
library(ade4)
library(adegenet)
library(ape)
library(pegas)
library(doParallel)
library(foreach)
library(dplyr)

# read the utility functions
source("../../../../../study_pipeline/workflow/scripts/utilities.R")


# sliding window Tajima D function
tajima_d_calculation <- function(aln_rec_chr, outgroup, pop_metadata, coords_file, output_table) {

    # reaed the meta data
    pop_meta <- population_host_metadata(pop_metadata, host = FALSE, humans = TRUE)

    # read the alignment for the sliding window Tajima D
    aln <- read_aln(aln_rec_chr, outgroup, sw = TRUE, tree = FALSE)

    # read the file with the coordinates for each gene
    annot <- read.csv(coords_file, sep='\t', header=TRUE)

    # define the dataframe where were storing the results
    tajima_d_df <- data.frame(matrix(ncol = 3, nrow = 0))

    # start the SW
    for (i in annot$ID) {
        print(paste("Analysing gene", i, sep=' '))
        gene <- i
        start <- annot %>% filter(ID == i) %>% getElement('start')
        end <-annot %>% filter(ID == i) %>% getElement('end')
        chunk <- aln[,start:end]

        for (pop in unique(pop_meta$Trait)) {
            d <- tajima.test(chunk[pop_meta[pop_meta$Trait == pop,'sample'],])$D
            tajima_d_df <- rbind.data.frame(tajima_d_df, c(pop, round(d, 2), gene))
        }
        names(tajima_d_df) <- c('Population', 'Tajima D', 'Bin')
        # store an Fst matrix if successful, otherwise store NULL
     }

    write.table(tajima_d_df, file=output_table, row.names=FALSE, quote=FALSE, sep='\t')
}

# run the function
tajima_d_calculation(aln_rec_chr, outgroup, pop_metadata, coords_file, output_table)