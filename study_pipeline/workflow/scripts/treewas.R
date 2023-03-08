# Title     : treeWAS
# Objective : Run TreeWAS
# Created by: Evangelos A. Dimopoulos"
# Created on: 08/04/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
outgroup <- args[2]
aln_file <- args[3]
pop_metadata <- args[4]
plot_file <- args[5]
out_table <- args[6]

# load libs
library(treeWAS)

# read the utility functions
source("scripts/utilities.R")

# treeWAS function
treewas <- function(tree_file, outgroup, aln_file, pop_metadata, plot_file, out_table) {

    # read tree
    tree_r_no_out <- read_tree(tree_file, outgroup)

    # read aln data from file:
    genind_obj <- read_aln(aln_file, tree_r_no_out)
    mat <- genind_obj@tab

    # read file with population metadata
    pop_meta <- population_host_metadata(pop_metadata)

    # load traits
    pop_trait <- pop_meta[pop_meta$sample %in% tree_r_no_out$tip.label,]
    # pop_trait$Num_Trait <- ifelse(pop_trait$Trait == "Human", 1, 0)
    phen_num <- setNames(pop_trait$Num_Trait, pop_trait$sample)

    # run treeWas
    treewas_out <- treeWAS(mat, phen_num, tree = tree_r_no_out, seed = 12345, plot.tree = TRUE,
        plot.null.dist = TRUE, plot.dist = TRUE, plot.manhattan = TRUE, filename.plot = plot_file,
        phen.type = "categorical", p.value = 0.05, correct.prop = TRUE)

    # count number of significant SNPs
    terminal_nsnps <- if (is.vector(treewas_out$terminal$sig.snps)) {0} else {nrow(treewas_out$terminal$sig.snps)}
    simultaneous_nsnps <- if (is.vector(treewas_out$simultaneous$sig.snps)) {0
                          } else {
                          nrow(treewas_out$simultaneous$sig.snps)}
    subsequent_nsnps <- if (is.vector(treewas_out$subsequent$sig.snps)) {0} else {nrow(treewas_out$subsequent$sig.snps)}

    # put the treewas summary stats into vectors and then build a df with them
    res_terminal <- c("Terminal", terminal_nsnps, round(treewas_out$terminal$sig.thresh[[1]], digits = 3))
    res_simultaneous <- c("Simultaneous", simultaneous_nsnps,
    round(treewas_out$simultaneous$sig.thresh[[1]], digits = 3))
    res_subsequent <- c("Subsequent", subsequent_nsnps, round(treewas_out$subsequent$sig.thresh[[1]], digits = 3))

    treewas_res <- as.data.frame(rbind(res_terminal, res_simultaneous, res_subsequent))
    names(treewas_res) <- c("TreeWAS test", "Num of Sig SNPs", "Sig Threshold")

    # write output table
    write.table(treewas_res, file=out_table, row.names=FALSE, quote=FALSE, sep='\t')
}

# run the function
treewas(tree_file, outgroup, aln_file, pop_metadata, plot_file, out_table)