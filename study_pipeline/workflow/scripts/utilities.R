# Title     : utilities
# Objective : R utility functions
# Created by: Evangelos A. Dimopoulos"
# Created on: 10/04/2022


root_tree <- function(tree_file, outgroup) {
  # read tree
  tree <- read.tree(file = tree_file)
  tree_r <- root(tree, outgroup = outgroup, edgelabel = TRUE)
  return(tree_r)
}


read_tree <- function(tree_file, outgroup) {
  # read tree
  tree <- read.tree(file = tree_file)
  tree_r <- root(tree, outgroup = outgroup, edgelabel = TRUE)
  tree_r_no_out <- drop.tip(tree_r, c(outgroup), rooted = TRUE, collapse.singles = TRUE)
  cat("The number of tips in the tree is: ", length(tree_r_no_out$tip.label), "\n")
  return(tree_r_no_out)
}


read_aln <- function(aln_file, set_no_out, sw = FALSE, tree = TRUE) {
  # read data from file:
  dna <- read.dna(file = aln_file, format = "fasta", as.matrix = TRUE)
  dna <- if (tree) dna[set_no_out$tip.label,] else dna[!rownames(dna) %in% c(set_no_out),]
  cat("The number of samples in the alignment is: ", dim(dna)[1], "\n")
  cat("the total length of the alignment is: ", dim(dna)[2], "\n")

  # convert dna to genid and get the nucleotide matrix
  genid_obj <- DNAbin2genind(dna)

  # based on if I am doing a sliding window either return the DNAbin object or the Genind object
  aln_obj <- if (sw) dna else genid_obj
  return(aln_obj)
}


population_host_metadata <- function(pop_metadata, host = TRUE, humans = FALSE) {
  # read file with population metadata
  pop_meta <- read.csv(file = pop_metadata, sep = "\t", header = TRUE)

  if (host) {
    human <- c("Human")
    ruminant <- c("Beef", "Bovine", "Sheep", "Goat", "Ovine/Goat")
    avian <- c("Chicken", "Avian", "Crow", "Poultry", "Turkey")
    food <- c("Food", "Dairy")
    swine <- c("Pork", "Swine", "Pig")
    other_mammal <- c("Primate", "Rodent", "Deer", "Canine", "Feline", "Equine", "Marine Mammal", "Other Mammal")
    other <- c("Soil", "Water", "Water/River", "Soil/Dust", "ND/Other", "Laboratory", "Plant", "Animal-related",
               "Environment")

    pop_meta$Trait[pop_meta$Host %in% human] <- "Human"
    pop_meta$Trait[pop_meta$Host %in% ruminant] <- "Ruminant"
    pop_meta$Trait[pop_meta$Host %in% avian] <- "Avian"
    pop_meta$Trait[pop_meta$Host %in% food] <- "Food"
    pop_meta$Trait[pop_meta$Host %in% swine] <- "Swine"
    pop_meta$Trait[pop_meta$Host %in% other_mammal] <- "Other_mammal"
    pop_meta$Trait[pop_meta$Host %in% other] <- "Other"

    trait_num_code <- setNames(c("Human", "Ruminant", "Avian", "Food", "Swine", "Other_mammal", "Other"),
                               c(1, 2, 3, 4, 5, 6, 7))
  } else if (humans) {
    pop_meta$Trait <- "Rest_of_Context"
    pop_meta[(pop_meta$Host == "Human") & (pop_meta$Dataset == "SB27"), 'Trait'] <- 'SB27_Humans'
    pop_meta[(pop_meta$Host == "Human") & (pop_meta$Dataset == "Context") &
               (pop_meta$Country == "USA"), 'Trait'] <- 'Context_America_Humans'

    trait_num_code <- setNames(c("SB27_Humans", "Context_America_Humans", "Rest_of_Context"),
                               c(1, 2, 3))
  }

  pop_meta$Num_Trait <- names(trait_num_code)[match(pop_meta$Trait, trait_num_code)]

  return(pop_meta)
}


parse_eucd <- function(eucd_df, metadata_df, dist_type, cohort) {
  # subset the euclidean distribution dataframe to only keep SB27 (rows) vs Context samples (cols)
  meta <- population_host_metadata(metadata_df)

  if (cohort == 'sb27') {
    set_1 <- "SB27"
    set_2 <- "Context"
    meta_set1 <- meta[meta$Dataset == set_1,]
    meta_set2 <- meta[meta$Dataset == set_2,]
  } else if (cohort == 'host') {
    set_1 <- "Human"
    set_2 <- "Animal"
    meta_set1 <- meta[meta$Host == set_1,]
    meta_set2 <- meta[!meta$Trait %in% c('Human', 'Food', 'Other'),]
  } else if (cohort == 'humans') {
    set_1 <- "SB27"
    set_2 <- "Human_context"
    meta_set1 <- meta[meta$Dataset == set_1,]
    meta_set2 <- meta[(meta$Host == "Human") & (meta$Dataset == "Context"),]
  } else if (cohort == 'americans') {
    set_1 <- "SB27"
    set_2 <- "Human_american"
    meta_set1 <- meta[meta$Dataset == set_1,]
    meta_set2 <- meta[(meta$Host == "Human") &
                        (meta$Dataset == "Context") &
                        (meta$Country == "USA") &
                        (meta$continent == "North America"),]
  } else if (cohort == 'americans-time') {
    set_1 <- "SB27"
    set_2 <- "Human_american_2015_2016"
    meta_set1 <- meta[meta$Dataset == set_1,]
    meta_set2 <- meta[(meta$Host == "Human") &
                        (meta$Dataset == "Context") &
                        (meta$Country == "USA") &
                        (meta$continent == "North America") &
                        (meta$year %in% c(2015, 2016)),]
  } else if (cohort == 'america-all') {
    set_1 <- "SB27"
    set_2 <- "American_context"
    meta_set1 <- meta[meta$Dataset == set_1,]
    meta_set2 <- meta[(meta$Dataset == "Context") & (meta$Country == "USA") & (meta$continent == "North America"),]
  } else if (cohort == 'narms-2019') {
    set_1 <- "Cali_NARMS_2019"
    set_2 <- "Rest_USA_NARMS_2019"
    meta_set1 <- meta[(meta$Dataset == "Context") &
                        (meta$year == 2019) &
                        (meta$Country == "USA") &
                        (meta$State == "USA:CA"),]
    meta_set2 <- meta[(meta$Dataset == "Context") &
                        (meta$year == 2019) &
                        (meta$Country == "USA") &
                        !(meta$State %in% c("USA:CA", "")),]
  } else if (cohort == 'sb27-cali-narms') {
    set_1 <- "SB27"
    set_2 <- "Cali_NARMS_2019"
    meta_set1 <- meta[meta$Dataset == set_1,]
    meta_set2 <- meta[(meta$Dataset == "Context") &
                        (meta$year == 2019) &
                        (meta$Country == "USA") &
                        (meta$State == "USA:CA"),]
  } else if (cohort == 'sb27-rest-narms') {
    set_1 <- "SB27"
    set_2 <- "Rest_USA_NARMS_2019"
    meta_set1 <- meta[meta$Dataset == set_1,]
    meta_set2 <- meta[(meta$Dataset == "Context") &
                        (meta$year == 2019) &
                        (meta$Country == "USA") &
                        !(meta$State %in% c("USA:CA", "")),]
  } else {
    stop("Options for this sample cohort haven't been coded yet. Either use sb27 or source.")
  }

  eucd_subset <- eucd_df[rownames(eucd_df) %in% meta_set1$sample, colnames(eucd_df) %in% meta_set2$sample]
  if (dim(eucd_subset)[1] == 0) {
    eucd_subset <- eucd_df[rownames(eucd_df) %in% meta_set2$sample, colnames(eucd_df) %in% meta_set2$sample]
    dist_type <- paste(dist_type, paste('No', set_1, sep = ' '), sep = '-')
  }
  if (dim(meta_set2)[1] == 0) {
    eucd_subset <- eucd_df[rownames(eucd_df) %in% meta_set1$sample, colnames(eucd_df) %in% meta_set1$sample]
    dist_type <- paste(dist_type, paste('No', set_2, sep = ' '), sep = '-')
  }
  eucd_subset <- rownames_to_column(eucd_subset, 'sample')

  # melt the dataframe
  eucd_melted <- melt(eucd_subset, id.vars = "sample")
  eucd_melted$dist <- dist_type
  return(eucd_melted)
}


big_exponent_limits <- function(x) {
  # set too big or too small numbers to 1.e308 or -1.e308
  x[x > 1e+308] <- 1e+308
  x[x < 1e-308] <- 1e-308
  return(x)
}