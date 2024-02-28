# Title     : pop gen stats
# Objective : Run Nei's (1982) and Weir and Cockerham's Fst, Tajima's D, Nei's Dxy, and broad sense heritability H^2
# Created by: Evangelos A. Dimopoulos"
# Created on: 08/04/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
outgroup <- args[1]
aln_file <- args[2]
pop_metadata <- args[3]
metric <- args[4]
output_table <- args[5]

# load libs
library(ade4)
library(adegenet)
library(ape)
library(hierfstat)
library(pegas)
library(reshape2)
library(poppr)

# read the utility functions
source("scripts/utilities.R")

# define the parameters for the sliding window scans
window_size <- 1000
overlap <- 900
step <- window_size - overlap

# sliding window Tajima D function
sw_tajima_d_scan <- function(aln_rec_chr, outgroup, pop_meta) {

    # read the alignment for the sliding window Tajima D
    aln <- read_aln(aln_rec_chr, outgroup, sw = TRUE, tree = FALSE)

    # starting sequence
    starts <- if (dim(aln)[2] > window_size) seq(1, dim(aln)[2], by = step) else seq(1:1)
    n <- length(starts)
    tajima_d_df <- data.frame(matrix(ncol = 3, nrow = 0))

    # start the SW
    for (i in 1:n) {
        start <- if (length(starts) > 1) (starts[i]) else 1
        end <- if ((start + step - 1) < dim(aln)[2]) (start + step - 1) else dim(aln)[2]
        chunk <- aln[,start:end]
        genome_bin <- if (length(starts) > 1) (starts[i]) else 1

        for (pop in unique(pop_meta$Trait)) {
            d <- tajima.test(chunk[pop_meta[pop_meta$Trait == pop,'sample'],])$D
            tajima_d_df <- rbind.data.frame(tajima_d_df, c(pop, round(d, 2), genome_bin))
        }
        names(tajima_d_df) <- c('Population', 'Tajima D', 'Bin')
        # store an Fst matrix if successful, otherwise store NULL
     }

    return(tajima_d_df)
}

# sliding window WC Fst function
sw_fst_scan <- function(aln_rec_chr, outgroup, pop_meta) {

    # read the alignment for the sliding window Fst
    aln <- read_aln(aln_rec_chr, outgroup, sw = TRUE, tree = FALSE)

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
        tryCatch({chunk_res <- hierfstat_calculate(outgroup, chunk, pop_meta, "swfst")},
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
hierfstat_calculate <- function(outgroup, dna_obj, pop_meta, metric) {

    # load traits IN THE CORRECT ORDER - same order as the dna_obj
    pop_trait <- pop_meta[pop_meta$sample %in% indNames(dna_obj[indNames(dna_obj) != outgroup]),]
    pop_trait <- pop_trait[match(indNames(dna_obj), pop_trait$sample),]
    strata(dna_obj) <- pop_trait[, c("Trait", "Host", "Source")]
    setPop(dna_obj) <- ~ Trait

    # convert genind object to hierfstat input
    hierfstat_input <- genind2hierfstat(dna_obj)

    # if sw or pairwise fst only calculate WC's pairwise Fst otherwise get the right metric
    if (metric == 'dxy') {
        # calculate Nei's Ds - the 1984 formulation also known as Da
        nei_dist <- genet.dist(hierfstat_input, diploid=FALSE, method="Da")

        # store the at the right slot
        stats <- list(NULL, NULL, nei_dist)
    } else if (metric == 'fst') {
        # calculate Nei's Fst - based on genotype frequencies
        nei_fstat <- basic.stats(hierfstat_input, diploid=FALSE,digits=2)$overall

        # calculate WC's Fst - based on ANOVA and DNA sequence data
        wc_fstat <- wc(hierfstat_input, diploid=FALSE)

        # store the at the right slot
        stats <- list(c(nei_fstat['Fst'], wc_fstat$FST), NULL, NULL)
    } else {
        # calculate WC's pairwise Fst - we prefer this to Nei's as the latter overestimates structure
        wc_pair_fstat <- pairwise.WCfst(hierfstat_input, diploid=FALSE)

        # store the at the right slot
        stats <- list(NULL, wc_pair_fstat, NULL)
    }

    return(stats)
}

# poppr heritability function
heritability_calculate <- function(outgroup, aln_file, pop_meta) {

    # read aln data from file:
    genind_obj <- read_aln(aln_file, outgroup, tree = FALSE)

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

# function to calculate Fsts/Dxy/TajiD/H^2 before and after masking recombination
pop_gen_stats  <- function(outgroup, aln_file, pop_metadata, metric, output_table) {

    # read file with population metadata
    if (metric == 'swfst-human' | metric == 'tajima-human') {
        # the trait is human focused
        pop_meta <- population_host_metadata(pop_metadata, host = FALSE, humans = TRUE)
    } else {
        pop_meta <- population_host_metadata(pop_metadata)
    }

    # run the actual pop gen stats calculation
    if (metric == 'swfst' | metric == 'swfst-human') {
        # run sliding window fst scan across the genome
        res_dist <- sw_fst_scan(aln_file, outgroup, pop_meta)
    } else if (metric == 'tajima' | metric == 'tajima-human') {
        res_dist <- sw_tajima_d_scan(aln_file, outgroup, pop_meta)
    } else if (metric == 'heritability') {
        res_dist <- heritability_calculate(outgroup, aln_file, pop_meta)
    } else {
        # distance results for a given alignment - order (Nei and WC total Fst, WC pair Fst, Nei' Dxy
        res_dist <- hierfstat_calculate(outgroup, read_aln(aln_file, outgroup, tree = FALSE), pop_meta, metric)
    }

    # depending on the stat create the appropriate output
    if (metric == 'dxy') {
        # write Nei's Ds output table
        write.table(as.matrix(res_dist[[3]]), file=output_table, quote=FALSE, row.names=TRUE, sep="\t",
            append=TRUE)
    } else if (metric == 'tajima' | metric == 'tajima-human') {
        # write the SW tajima D scan results
        write.table(res_dist, file=output_table, row.names=FALSE, quote=FALSE, sep='\t')
    } else if (metric == 'fst') {
        # check if the alignment has recombinant sites masked
        rec_state <- if (grepl('nrec', aln_file)) "No Recombination" else "Recombination"
        fst_res <- as.data.frame(t(c(rec_state, res_dist[[1]])))
        names(fst_res) <- c("Rec state", "Nei's Fst", "WC's Fst")

        # write Fst summary output table
        write.table(fst_res, file=output_table, row.names=FALSE, quote=FALSE, sep='\t')
    } else if (metric == 'pairfst') {
        # write the pairwise WC Fst
        write.table(as.matrix(res_dist[[2]]), file=output_table, quote=FALSE, row.names=TRUE, sep="\t",
        append=TRUE)
    } else if (metric == 'swfst' | metric == 'swfst-human') {
        write.table(res_dist, file=output_table, row.names=FALSE, quote=FALSE, sep='\t')
    } else if (metric == 'heritability') {
        # check if the alignment has recombinant sites masked
        rec_state <- if (grepl('nrec', aln_file)) "No Recombination" else "Recombination"
        h_broad_res <- as.data.frame(t(c(rec_state, res_dist)))
        names(h_broad_res) <- c("Rec state", "Heritability (broad sense)")

        # write output table
        write.table(h_broad_res, file=out_table, row.names=FALSE, quote=FALSE, sep='\t')
    }

}

# run the install function
pop_gen_stats(outgroup, aln_file, pop_metadata, metric, output_table)