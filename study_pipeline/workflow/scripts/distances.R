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
aln_file_nrec_chr <- args[5]
pop_metadata <- args[6]
sum_out_table <- args[7]
pair_wc_fst_out_table <- args[8]
pair_nei_dist_out_table <- args[9]
fst_sw_out_table <- args[10]

# load libs
library(ade4)
library(adegenet)
library(ape)
library(hierfstat)
library(pegas)
library(reshape2)

# read the utility functions
source("scripts/utilities.R")

# sliding window function

sw_fst_scan <- function(aln_file_nrec_chr, tree_obj, pop_meta) {

    # define the parameters for the sliding window Fst
    aln <- read_aln(aln_file_nrec_chr, tree_obj, sw = TRUE)
    window_size <- 10000
    overlap <- 9000
    step <- window_size - overlap

    # starting sequence
    starts <- if (dim(aln)[2] > window_size) seq(1, dim(aln)[2], by = step) else seq(1:1)
    n <- length(starts)
    genome_bins <- list()

    # start the SW
    for (i in 1:n) {
        start <- if (length(starts) > 1) (starts[i]) else 1
        end <- if ((start + step - 1) < dim(aln)[2]) (start + step - 1) else dim(aln)[2]
        chunk <- DNAbin2genind(aln[,start:end])
        chunk_res <- list(NULL, NULL, NULL)
        # if the segment has either samples where all the positions are NA or they are monomorphic then this will fail
        tryCatch({chunk_res <- hierfstat_calculate(tree_obj, chunk, pop_meta, sw = TRUE)},
            error= function(e) {cat("ERROR :",conditionMessage(e), "\n",
            "Fst cannot be calculated for that region. There is either too much missing data or no polymorphic positions.",
            "\n")})

        # store an Fst matrix if successful, otherwise store NULL
        if (is.null(chunk_res[[2]])) {
        genome_bins[i] <- list(NULL)
        } else {
        genome_bins[[i]] <- chunk_res[[2]]
        }
     }

    # store the results
    sw_fst <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Host_1", "Host_2", "Fst", "Bin"))
    if (length(genome_bins) > 0) {
    for (i in 1:n) {
        # if the element of the list is NULL it means that this bin has failed to produce an Fst
        if (is.null(genome_bins[[i]])) next

        # convert the matrix format to long dataframe format
        matrix_input <- genome_bins[[i]]
        matrix_input[is.na(matrix_input)] <- 0
        matrix_input[upper.tri(matrix_input)] <- NA
        fst_df_melt <- melt(matrix_input, na.rm = TRUE)
        fst_df_melt$bin <- if (length(starts) > 1) (starts[i]) else 1
        sw_fst <- rbind.data.frame(sw_fst, fst_df_melt)
    }
    names(sw_fst) <- c("Host_1", "Host_2", "Fst", "Bin")

    return(sw_fst)
    } else {
    return(sw_fst)
    }
}

# hierfstat function
hierfstat_calculate <- function(tree_obj, dna_obj, pop_meta, sw = FALSE) {

    # load traits
    pop_trait <- pop_meta[pop_meta$sample %in% tree_obj$tip.label,]
    phen_cat <- setNames(pop_trait$Trait, pop_trait$sample)

    # convert genind object to hierfstat input
    hierfstat_input <- genind2hierfstat(dna_obj, pop=phen_cat)

    # if not sliding windows calculate all the stats, otherwise do not waste time
    if (!sw) {
        # calculate Nei's Fst - based on genotype frequencies
        nei_fstat <- basic.stats(hierfstat_input, diploid=FALSE,digits=2)$overall

        # calculate WC's Fst - based on ANOVA and DNA sequence data
        wc_fstat <- wc(hierfstat_input, diploid=FALSE)

        # calculate Nei's Ds - the 1984 formulation also known as Da
        nei_dist <- genet.dist(hierfstat_input, diploid=FALSE, method="Da")
    }

    # calculate WC's pairwise Fst - we prefer this to Nei's as the latter overestimates structure
    wc_pair_fstat <- pairwise.WCfst(hierfstat_input, diploid=FALSE)

    # chose what list to return based on whether we do sliding windows or not
    if (sw) {
    stats <- list(NULL, wc_pair_fstat, NULL)
    } else {
    stats <- list(c(nei_fstat['Fst'], wc_fstat$FST), wc_pair_fstat, nei_dist)
    }

    return(stats)
}

# function to calculate Fsts before and after masking recombination
distances  <- function(tree_file, outgroup, aln_file_rec, aln_file_nrec, aln_file_nrec_chr, pop_metadata,
    sum_out_table, pair_wc_fst_out_table, pair_nei_dist_out_table, fst_sw_out_table) {

    # read tree
    tree_r_no_out <- read_tree(tree_file, outgroup)

    # read file with population metadata
    pop_meta <- population_host_metadata(pop_metadata)

    # get the Fsts for the alignments before we remove the recombination
    rec_dist <- hierfstat_calculate(tree_r_no_out, read_aln(aln_file_rec, tree_r_no_out), pop_meta)

    # get the Fsts for the alignments after we remove the recombination
    nrec_dist <- hierfstat_calculate(tree_r_no_out, read_aln(aln_file_nrec, tree_r_no_out), pop_meta)

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

    # run sliding window fst scan across the genome
    sw_fst <- sw_fst_scan(aln_file_nrec_chr, tree_r_no_out, pop_meta)

    write.table(sw_fst, file=fst_sw_out_table, row.names=FALSE, quote=FALSE, sep='\t')

}

# run the install function
distances(tree_file, outgroup, aln_file_rec, aln_file_nrec, aln_file_nrec_chr, pop_metadata,
    sum_out_table, pair_wc_fst_out_table, pair_nei_dist_out_table, fst_sw_out_table)