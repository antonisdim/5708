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
suppressMessages(library(rstatix))

# read the utility functions
source("../../../study_pipeline/workflow/scripts/utilities.R")

PREVALENCE <- 0.90

gene_carriage_comp <- function(pop_metadata, abricate_report, table_out_name, sb27, lineage) {

  # check this sb27 argument is correct
  if (!(sb27 %in% c('all', 'sb27', 'sb27_strict', 'sb27_inv'))) stop("It is either all or sb27")

  # read the pop metadata file
  pop_meta <- population_host_metadata(pop_metadata)

  # read the abridate report file
  gene_df <- read.csv(abricate_report, sep = '\t', header = TRUE)

  # substitute '.' with NA
  gene_df[gene_df == "."] <- NA

  # take away the core and soft core genes
  gene_df_melt <- melt(gene_df, na.rm = FALSE, value.name = "gene_id_percent",
                     variable.name = "gene", id = c("X.FILE", "NUM_FOUND"))

  # calculate the prevalnce of the amr genes in the dataset
  gene_prevalnce <- gene_df_melt %>%
    drop_na(gene_id_percent) %>%
    group_by(gene) %>%
    summarise(n=n()) %>%
    mutate(prevalence = n/n_distinct(gene_df_melt$X.FILE))

  gene_counts <- gene_df_melt %>%
    filter(gene %in% gene_prevalnce[gene_prevalnce$prevalence <= PREVALENCE,]$gene) %>%
    mutate(presence = ifelse(is.na(gene_id_percent), 0, 1)) %>%
    group_by(X.FILE) %>%
    summarise(real_found_genes = sum(presence))

  # join the two dataframes
  pop_meta <- pop_meta %>% right_join(gene_counts, by = c("sample" = "X.FILE"))

  # filter based on whether we are keeping only the SB27 samples or not and which lineage
  if (sb27 == 'sb27') pop_meta <- pop_meta[pop_meta$Dataset == 'SB27',]

  # strict mode - only sick humans from SB27 and animals that were both raised and purcahsed in California
  if (sb27 == 'sb27_strict') pop_meta <- pop_meta %>% filter(Dataset == 'SB27' & Source == 'Human' |
                                                               (Source == 'Animal' & CA == 1))

  # inverse means keep only animals bought in California that were imported
  if (sb27 == 'sb27_inv') pop_meta <- pop_meta %>% filter(Dataset == 'SB27' & Source == 'Animal' & CA == 0)

  # select a lineage
  if (as.numeric(lineage) != 1000) pop_meta <- pop_meta[pop_meta$cluster == as.numeric(lineage),]

  cat(paste("We are analysing", nrow(pop_meta), "isolates\n", sep = " "))

  # get a list of dates and keep only the ones after 2016
  pop_meta <- pop_meta[pop_meta$Collection_Year >= 2016, ]

  # store the results
  pvalues_df <- data.frame(matrix(ncol = 8, nrow = 0))

  # loop through the sources
  for (source in c('all', 'Human', 'Animal')) {

    # get the right subset of the dataset
    if (source != 'all') pop_test <- pop_meta[pop_meta$Source == source,] else pop_test <- pop_meta

    # if no samples skip to the next dataset source
    if (dim(pop_test)[1] == 0) next

    # compare the distributions of the amr/vf genes carried by isolates in every year
    # pairwise t-tests, not paired, bonferroni correction, alternative hypothesis is 'less'
    ttest_res <- pairwise_t_test(pop_test, real_found_genes ~ Collection_Year, p.adjust.method='bonferroni',
                                 alternative='less', paired = FALSE)

    ttest_res_sum <- ttest_res %>% select('group1', 'group2', 'n1', 'n2', 'p.adj', 'p.adj.signif')

    # add the lineage and the source labels
    ttest_res_sum$lineage <- as.numeric(lineage)
    ttest_res_sum$source <- source

    # append to the dataframe
    pvalues_df <- rbind.data.frame(pvalues_df, ttest_res_sum)

    # report results
    cat(paste('Comparing if the carriage of genes differs among yearly cohorts for', source, 'and cluster',
              lineage, "\n", sep = " "))

  }

  # rename the columns to something meaningful
  names(pvalues_df) <- c('year_1', 'year_2', 'n_1', 'n_2', 'p_adj', 'p_adj_signif', 'lineage', 'source')

  # write the table
  write.table(pvalues_df, file = table_out_name, quote = FALSE, row.names = FALSE, sep = "\t")
}

# run the function
gene_carriage_comp(pop_metadata, abricate_report, table_out_name, sb27, lineage)
