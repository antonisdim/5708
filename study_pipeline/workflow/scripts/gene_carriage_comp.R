# Title     : gene_carriage_comparison
# Objective : Compare the number of amr or vf genes that are carried by isolates each year
# Created by: Evangelos A. Dimopoulos"
# Created on: 16/11/2023

# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
pop_metadata <- args[1]
abricate_report <- args[2]
table_out_name <- args[3]
sb27 <- args[4]
lineage <- args[5]

# load libs
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(harmonicmeanp))
suppressMessages(library(ggplot2))

# read the utility functions
source("../../../study_pipeline/workflow/scripts/utilities.R")


gene_carriage_comp <- function(pop_metadata, abricate_report, table_out_name, sb27, lineage) {

  # check this sb27 argument is correct
  if (!(sb27 %in% c('all', 'sb27'))) stop("It is either all or sb27")

  # read the pop metadata file
  pop_meta <- population_host_metadata(pop_metadata)

  # read the abridate report file
  gene_counts <- read.csv(abricate_report, sep = '\t', header = TRUE)
  gene_counts <- gene_counts[, c('X.FILE', 'NUM_FOUND')]

  # join the two dataframes
  pop_meta <- pop_meta %>% right_join(gene_counts, by = c("sample" = "X.FILE"))

  # filter based on whether we are keeping only the SB27 samples or not and which lineage
  if (sb27 == 'sb27') pop_meta <- pop_meta[pop_meta$Dataset == 'SB27',]
  if (as.numeric(lineage) != 1000) pop_meta <- pop_meta[pop_meta$cluster == as.numeric(lineage),]
  cat(paste("We are analysing", nrow(pop_meta), "isolates", sep = " "))

  # get a list of dates and keep only the ones after 2016
  date_list <- sort(unique(pop_meta$Collection_Year))
  date_list <- date_list[date_list >= 2016]

  # loop through the pairs and use a counter to see when you've reached 2022 - that counter == the index of the maximum element
  max_date <- max(date_list)

  # store the results
  pvalues_df <- data.frame(source = character(), before = double(), after = double(),
                           Welch_pvalue = double(), Wilcox_pvalue = double(), lineage = double())

  for (year in 1:length(date_list)) {
    # we are effectively looping in pairs - so we are stopping at the second to last year
    if (date_list[year] < max_date) {

      # define before and after years
      before <- date_list[year]
      after <- date_list[year + 1]

      # subset the dataframe accrodingly
      for (source in c('all', 'Human', 'Animal')) {
        if (source == 'all') {
          isolates_before <- pop_meta[pop_meta$Collection_Year == before,]
          isolates_after <- pop_meta[pop_meta$Collection_Year == after,]
        } else {
          isolates_before <- pop_meta[(pop_meta$Collection_Year == before) & (pop_meta$Source == source),]
          isolates_after <- pop_meta[(pop_meta$Collection_Year == after) & (pop_meta$Source == source),]
        }

        # compare the distributions of the amr/vf genes carried by isolates in every year - they need to allow for vairable sample size
        # welch's t-test and
        # premeptive in case the first iteration raises an error
        welch <- list()
        tryCatch({ welch <- t.test(isolates_before$NUM_FOUND, isolates_after$NUM_FOUND) },
                 error = function(e) { cat("ERROR :", conditionMessage(e), "\n", "not enough data\n")
                   welch$p.value <<- NaN })
        # Mann-Whitney-Wilcoxon's test
        wilcox <- list()
        tryCatch({ wilcox <- wilcox.test(isolates_before$NUM_FOUND, isolates_after$NUM_FOUND) },
                 error = function(e) { cat("ERROR :", conditionMessage(e), "\n", "not enough data\n")
                   wilcox$p.value <<- NaN })

        # append to the p-values dataframe
        test_vector <- c(source, before, after, welch$p.value, wilcox$p.value, lineage)
        pvalues_df[nrow(pvalues_df) + 1,] <- test_vector

        # report the results
        cat(paste('Comparing if the carriage of genes differs between', before, 'and', after,
                  'for source', source, '. Pvalue for the Welch t-test is', welch$p.value, 'and for the Wolcoxon test is',
                  wilcox$p.value, "\n", sep = " "))
      }
    }
  }

  # write the table
  write.table(pvalues_df, file = table_out_name, quote = FALSE, row.names = FALSE, sep = "\t")
}

# run the function
gene_carriage_comp(pop_metadata, abricate_report, table_out_name, sb27, lineage)
