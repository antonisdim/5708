# Title     : utilities
# Objective : R utility functions
# Created by: Evangelos A. Dimopoulos"
# Created on: 10/04/2022


read_tree <- function(tree_file, outgroup) {
    # read tree
    tree <- read.tree(file = tree_file)
    tree_r <- root(tree, outgroup = outgroup, edgelabel = TRUE)
    tree_r_no_out <- drop.tip(tree_r, c(outgroup), rooted = TRUE, collapse.singles = TRUE)
    return(tree_r_no_out)
}


read_aln <- function(aln_file, tree_r_no_out) {
    # read data from file:
    dna <- read.dna(file = aln_file, format = "fasta", as.matrix = TRUE)
    dna <- dna[tree_r_no_out$tip.label,]

    # convert dna to genid and get the nucleotide matrix
    genid_obj <- DNAbin2genind(dna)
    return(genid_obj)
}


population_host_metadata <- function(pop_metadata) {
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

    return(pop_meta)
}