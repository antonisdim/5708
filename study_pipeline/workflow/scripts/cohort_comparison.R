# Title     : cohort_comparison
# Objective : Compare the euclidean distances of cohorts of samples
# Created by: Evangelos A. Dimopoulos"
# Created on: 02/01/2023

# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
outgroup <- args[2]
aln_file <- args[3]
pop_metadata <- args[4]
hist_out_figure <- args[5]
pval_out_table <- args[6]
cohort <- args[7]

# load libs
library(adegenet)
library(ape)
library(tibble)
library(reshape2)
library(ggplot2)
library(harmonicmeanp)

# read the utility functions
source("scripts/utilities.R")

# palette
titipounamu <- c("#3E4331", "#AD6B17", "#66743B", "#D0C471", "#CCB62F", "#BAC4C2")

cohort_comparison <- function(tree_file, outgroup, aln_file, pop_metadata, hist_out_figure, pval_out_table, cohort) {

  # read tree
  tree_r_no_out <- read_tree(tree_file, outgroup)

  # read aln data from file:
  genind_obj <- read_aln(aln_file, tree_r_no_out)

  # get the euclidean distances
  euclidean_dist <- dist(genind_obj, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)

  # turn the dist object into a melted dataframe
  eucd_df <- as.data.frame(as.matrix(euclidean_dist, labels = TRUE))
  eucd_obs <- parse_eucd(eucd_df, pop_metadata, "Observed distribution", cohort)

  # define an empty dataset to store p-values and randomised euclidean distance datasets
  pvalues_df <- data.frame(iter = character(), Welch_pvalue = double(), Wilcox_pvalue = double())
  randomised_distances_list <- list()

  # in a for loop this
  for (iter in 1:100) {
    eucd_iter <- eucd_df
    colnames(eucd_iter) <- sample(colnames(eucd_iter))
    eucd_rand <- parse_eucd(eucd_iter, pop_metadata, "Random sample", cohort)

    # calculate welch's t-test and Mann-Whitney-Wilcoxon's test
    welch <- t.test(eucd_obs$value, eucd_rand$value)
    wilcox <- wilcox.test(eucd_obs$value, eucd_rand$value)

    test_vector <- c(iter, welch$p.value, wilcox$p.value)
    pvalues_df[nrow(pvalues_df) + 1,] <- test_vector

    # for plots
    compare <- rbind(eucd_rand, eucd_obs)
    compare$iter <- iter

    randomised_distances_list[[iter]] <- compare
  }

  # construct dataframe of iterations to plot
  distances_comparison <- do.call(rbind, randomised_distances_list)

  eucd_comparison_histogram <- ggplot(data = distances_comparison, aes(x = value, fill = dist)) +
    geom_histogram(alpha = 0.9) +
    theme_minimal() +
    xlab("Euclidean distance") +
    scale_fill_manual(values = titipounamu) +
    theme(plot.title = element_text(size = 12, hjust = 0.5)) +
    facet_wrap(~iter)

  ggsave(hist_out_figure, eucd_comparison_histogram, width = 12, height = 12)

  #calculate the harmonic mean of the p-values from each iteration
  pval_total_vector <- c("total", p.hmp(big_exponent_limits(pvalues_df$Welch_pvalue), w = NULL, L = 100),
                         p.hmp(big_exponent_limits(pvalues_df$Wilcox_pvalue), w = NULL, L = 100))
  pvalues_df[nrow(pvalues_df) + 1,] <- pval_total_vector

  write.table(pvalues_df, file = pval_out_table, quote = FALSE, row.names = FALSE, sep = "\t")
}

# run the function
cohort_comparison(tree_file, outgroup, aln_file, pop_metadata, hist_out_figure, pval_out_table, cohort)

