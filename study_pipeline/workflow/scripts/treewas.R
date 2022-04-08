# Title     : r_github
# Objective : Install r packages that are not on conda
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

#load libs
library(treeWAS)

#treeWAS function
treewas <- function(tree_file, outgroup, aln_file, pop_metadata, plot_file, out_table) {

    # read tree
    tree <- read.tree(file = tree_file)
    tree_r <- root(tree, outgroup = outgroup, edgelabel = TRUE)
    tree_r_no_out <- drop.tip(tree_r, c(outgroup), rooted = TRUE, collapse.singles = TRUE)

    # read data from file:
    dna <- read.dna(file = aln_file, format = "fasta", as.matrix = TRUE)
    dna <- dna[tree_r_no_out$tip.label,]

    # convert dna to genid and get the nucleotide matrix
    genid_obj <- DNAbin2genind(dna)
    mat <- genid_obj@tab

    # read file with population metadata
    pop_meta <- read.csv(file = pop_metadata, sep = "\t")

    human <-c("Human")
    ruminant <- c("Beef", "Bovine", "Ovine/Goat")
    avian <- c("Chicken", "Avian", "Poultry", "Turkey")
    food <- c("Food", "Dairy")
    swine <- c("Pork", "Swine")
    other_mammal <- c("Primate", "Rodent", "Deer", "Canine")
    other <- c("Water/River", "Soil/Dust", "ND/Other", "Laboratory", "Plant", "Animal-related")

    pop_meta$Trait[pop_meta$Host %in% human] <- "Human"
    pop_meta$Trait[pop_meta$Host %in% ruminant] <- "Ruminant"
    pop_meta$Trait[pop_meta$Host %in% avian] <- "Avian"
    pop_meta$Trait[pop_meta$Host %in% food] <- "Food"
    pop_meta$Trait[pop_meta$Host %in% swine] <- "Swine"
    pop_meta$Trait[pop_meta$Host %in% other_mammal] <- "Other_mammal"
    pop_meta$Trait[pop_meta$Host %in% other] <- "Other"

    trait_num_code <- setNames(c("Human", "Ruminant", "Avian", "Food", "Swine", "Other_mammal", "Other"),
        c(1, 2, 3, 4, 5, 6, 7))

    pop_meta$Num_Trait <- names(trait_num_code)[match(pop_meta$Trait, trait_num_code)]

    # load traits
    pop_trait <- pop_meta[pop_meta$sample %in% tree_r_no_out$tip.label,]
    phen_num <- setNames(pop_trait$Num_Trait, pop_trait$sample)

    # run treeWas
    treewas_out <- treeWAS(mat, phen_num, tree = tree_r_no_out, seed = 12345, plot.tree = TRUE,
        plot.null.dist = TRUE, plot.dist = TRUE, plot.manhattan = TRUE, filename.plot = plot_file)

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

# run the install function
treewas(tree_file, outgroup, aln_file, pop_metadata, plot_file, out_table)