# Title     : cohort_comparison
# Objective : Compare the euclidean distances of cohorts of samples
# Created by: Evangelos A. Dimopoulos"
# Created on: 02/01/2023

# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
outgroup <- args[1]
in_file <- args[2]
pop_metadata <- args[3]
cohort <- args[4]
in_type <- args[5]

# load libs
suppressMessages(library(adegenet))
suppressMessages(library(ape))
suppressMessages(library(tibble))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))

# read the utility functions
source("../../../study_pipeline/workflow/scripts/utilities.R")

# palette
titipounamu <- c("#3E4331", "#AD6B17", "#66743B", "#D0C471", "#CCB62F", "#BAC4C2")

# number of permutation
n_perm <- 100

cohort_comparison <- function(outgroup, in_file, pop_metadata, cohort, in_type) {

  # read aln data from file:
  if (in_type == 'alignment') {
    genind_obj <- read_aln(in_file, outgroup, tree = FALSE)

    # get the euclidean distances
    euclidean_dist <- dist(genind_obj, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)

    # turn the dist object into a melted dataframe
    dist_df <- as.data.frame(as.matrix(euclidean_dist, labels = TRUE))

  } else {
     # read a distance matrix straight away
    dist_df <- read.csv(in_file, sep='\t', header=TRUE) %>% select(-X.)

    # redifine the row names
    rownames(dist_df) <- names(dist_df)
  }

  # get the dataframe for the observed distribution
  dist_obs <- parse_eucd(dist_df, pop_metadata, "Observed distribution", cohort, in_type)

  # get the observed D(a, b)
  observed_dist <- mean(as.matrix(dist_obs %>% select(-sample)))

  # define a list to store the permuted D(a, b)
  perm_dists <- numeric(n_perm)

  # in a for loop this
  for (iter in 1:100) {
    dist_iter <- dist_df
    colnames(dist_iter) <- sample(colnames(dist_iter))
    rownames(dist_iter) <- sample(rownames(dist_iter))
    dist_rand <- parse_eucd(dist_iter, pop_metadata, "Permuted distribution", cohort, in_type)

    perm_dists[iter] <- mean(as.matrix(dist_rand %>% select(-sample)))
    }

  # Calculate empirical p-value
  p_value_empirical <- sum(perm_dists >= observed_dist) / n_perm

  # use randtest to get a p-value
  rand_test_res <- as.randtest(sim = perm_dists, obs = observed_dist, alter = "greater")
  p_value_randtest <- rand_test_res$pvalue

  # print the results
  cat(paste("Cohort", "D(a, b) obs", "D(a, b) perm", "D(a, b) p-value", sep='\t'), "\n")
  cat(paste(cohort, observed_dist, mean(perm_dists), round(p_value_randtest, 2), sep='\t'), "\n")
}

# run the function
cohort_comparison(outgroup, in_file, pop_metadata, cohort, in_type)

