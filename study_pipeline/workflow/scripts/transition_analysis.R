# Title     : transition analysis
# Objective : Run host transition (switch) analysis
# Created by: Evangelos A. Dimopoulos"
# Created on: 23/06/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
snp_file <- args[1]
treefile <- args[2]
outgroup <- args[3]
pop_meta <- args[4]
all_out_table <- args[5]
most_common_table <- args[6]
state_sum_table <- args[7]
state_tr_table <- args[8]
clock <- as.numeric(args[9])
genome_size <- as.numeric(args[10])
time_interval <- as.numeric(args[11])
cluster <- as.numeric(args[12])

# load libs
library(igraph)
library(ape)
library(plyr)
library(dplyr)
library(tidygraph)
library(phytools)
library(tibble)

# read the utility functions
source("scripts/utilities.R")

# constants
CLOCK_STDEV <- 1.87e-7

snp_dist_net <-
  function(snp_dist,
           tree_obj,
           pop_meta,
           out_table,
           clock,
           genome_size,
           time_interval) {

    # read the snp distance file
    ecoli_pwsnps <- read.csv(snp_dist, sep = '\t', header = TRUE)

    # convert pairsnp output into dataframe
    ecoli_pwsnps_idx <- ecoli_pwsnps
    rownames(ecoli_pwsnps_idx) <- ecoli_pwsnps_idx$X.
    ecoli_pwsnps_idx <-
      ecoli_pwsnps_idx[, -which(names(ecoli_pwsnps_idx) %in% c("X", "X."))]

    pwsnps_decon <-
      data.frame(t(combn(names(ecoli_pwsnps_idx), 2)), dist = t(ecoli_pwsnps_idx)[lower.tri(ecoli_pwsnps_idx)])
    colnames(pwsnps_decon) <- c("Taxon1", "Taxon2", "dist")
    pwsnps_decon$Taxon1 <-
      gsub('# ', "", pwsnps_decon$Taxon1)

    # read the corresponding tree, create a network out of it and then convert it to df
    tree_no_out <- tree_obj

    tree_no_out$node.label <- NULL
    xg <- as.igraph(tree_no_out, directed = FALSE)
    comm <- cluster_fast_greedy(xg)
    V(xg)$community <- membership(comm)
    cat("Done with converting the tree to a network.\n")

    vcom <- data.frame(as_ids(V(xg)), V(xg)$community)
    colnames(vcom) <- c("sample", "network_id")

    network_phylo <- vcom[-grep("Node", vcom$sample),]

    # add host trait
    sample_meta_df <-
      population_host_metadata(pop_meta)

    # match sample and host
    pwsnps_decon$Taxon1_host <-
      sample_meta_df$Trait[match(pwsnps_decon$Taxon1, sample_meta_df$sample)]
    pwsnps_decon$Taxon2_host <-
      sample_meta_df$Trait[match(pwsnps_decon$Taxon2, sample_meta_df$sample)]

    # match sample and trial area
    pwsnps_decon$Taxon1_area <-
      sample_meta_df$Country[match(pwsnps_decon$Taxon1, sample_meta_df$sample)]
    pwsnps_decon$Taxon2_area <-
      sample_meta_df$Country[match(pwsnps_decon$Taxon2, sample_meta_df$sample)]

    # match sample and network_id
    pwsnps_decon$Taxon1_network <-
      network_phylo$network_id[match(pwsnps_decon$Taxon1, network_phylo$sample)]
    pwsnps_decon$Taxon2_network <-
      network_phylo$network_id[match(pwsnps_decon$Taxon2, network_phylo$sample)]

    # create new df with columns for area_match and host_match to help with separating and filtering data for plots
    pwsnps_decon_filtered <- pwsnps_decon %>%
      filter(Taxon1_area != "Unknown" &
               Taxon2_area != "Unknown") %>%
      mutate(area_match = ifelse(
        Taxon1_area == Taxon2_area,
        paste("same_trial_area"),
        paste("different_trial_area")
      )) %>%
      mutate(host_match = ifelse(
        Taxon1_host == Taxon2_host,
        paste("same_host"),
        paste("different_host")
      )) %>%
      mutate(network_match = ifelse(
        Taxon1_network == Taxon2_network,
        paste("Same Transmission Cluster"),
        paste("Different Transmission Cluster")
      )) %>%
      mutate(host_link = paste(Taxon1_host, Taxon2_host, sep = "-"))

    # calculate the SNP cutoff for the samples belonging to the same transmission group - diveded by 2 because of colaescence
    if (clock < 1e-05) {
      cat("Consider increasing the time interval if you don't see many transmission chains\n")
    }
    cat(paste("The clock rate is", clock, "\n", sep=" "))
    snp_cutoff <- round((clock * time_interval * genome_size) / 2)
    cat(paste("The snp cutoff is", snp_cutoff, "\n", sep=" "))

    # create new df filtering out samples that aren't in the same network and keep only the transmission pairs
    pwsnps_decon_filtered_clusters <-
      pwsnps_decon_filtered %>%
        filter(!is.na(Taxon1_host)) %>%
        filter(!is.na(Taxon2_host)) %>%
        filter(!is.na(network_match))
    # filter(!grepl("Other", host_link))

    pwsnps_decon_filtered_clusters_mini_pairs <-
      pwsnps_decon_filtered_clusters %>%
        filter(dist <= snp_cutoff) %>%
        filter(network_match == "Same Transmission Cluster")

    trans_summary <-
      pwsnps_decon_filtered_clusters_mini_pairs %>% count(host_link, sort = TRUE)

    write.table(
      trans_summary,
      file = out_table,
      row.names = FALSE,
      quote = FALSE,
      sep = '\t'
    )

  }

transition_analysis <-
  function(snp_dist,
           treefile,
           outgroup,
           pop_meta,
           all_out_table,
           most_common_table,
           state_sum_table,
           state_tr_table,
           clock,
           genome_size,
           time_interval,
           cluster) {


    tree_obj <- if (cluster == 1000) read_tree(treefile, outgroup) else read.nexus(treefile)
    snp_dist_net(snp_file, tree_obj, pop_meta, all_out_table, clock, genome_size, time_interval)
    cat("Done with transmission analysis for the whole tree.\n")

    sample_meta_df <- population_host_metadata(pop_meta)
    if (cluster != 1000) {
    sample_meta_df <- sample_meta_df[sample_meta_df$cluster == cluster,]
    }
    host_counts <- sample_meta_df %>%
      count(Trait, sort = TRUE) %>%
      arrange(n)
    subsample_num <- host_counts[3, c('n')]
    print.data.frame(host_counts)
    cat(paste("We are subsampling", subsample_num, "individuals from each category.\n", sep=" "))

    # hold the data to average over in these lists
    mean_state_time <- data.frame(matrix(ncol = nrow(sample_meta_df %>% count(Trait, sort = TRUE)) + 1, nrow = 0))
    total_state_changes <- c()
    total_tr_counts <- list()

    # create a directory to store the temporary host link count tables
    temp_dir <- paste(dirname(all_out_table), "/iter_tables_", sub('\\_pairsnp.tsv$', '', basename(snp_dist)), sep = "")

    dir.create(temp_dir)

    # start 100 iterations
    for (iter in 1:100) {
      out_table_iter <-
        paste(temp_dir,
              "/",
              sub('\\_pairsnp.tsv$', '', basename(snp_dist)),
              "_iter_",
              iter,
              ".tsv",
              sep = "")

      subsample_list <- list()
      cat(paste("Subsampling iteration", iter, "\n", sep=" "))

      # subsample the initital dataframe based with nsample equal to the second least frequent host type
      for (host in unique(sample_meta_df$Trait)) {
        host_count <- nrow(sample_meta_df[sample_meta_df$Trait == host,])

        # here is we need this otherwise an error is raised for nsample for the least frequent host type
        if (host_count < subsample_num) {
          nsample <- host_count
        } else {
          nsample <- subsample_num
        }

        # do the actual subsampling and storing to a list
        host_subsample <-
          sample_meta_df %>%
            filter(Trait == host) %>%
            sample_n(., nsample)
        subsample_list[[host]] <- host_subsample

      }

      # construct dataframe of the random set of samples
      subsample_df <- do.call(rbind, subsample_list)

      # drop the outgroup
      to_drop <-
        tree_obj$tip.label[!tree_obj$tip.label %in% subsample_df$sample]
      subset_tree <- tree_obj

      # prune the initial and keep only the samples that were randomly selected
      for (tip in to_drop) {
        subset_tree <- drop.tip(subset_tree, c(tip), rooted = TRUE, collapse.singles = TRUE)
      }

      # run the snp distance function
      snp_dist_net(snp_file, subset_tree, pop_meta, out_table_iter, clock, genome_size, time_interval)

      char <- unlist(as.list(deframe(subsample_df[, c("sample", "Trait")])))
      char <- char[names(char) != "SRR10270781"]
      print(subset_tree)
      # make the ACE simmap and summarise it
      stochastic_ace_tree <- make.simmap(subset_tree, char, model = "ER")
      ace_tree_sum <- describe.simmap(stochastic_ace_tree, plot = FALSE)

      # store the total state changes
      total_state_changes <- c(total_state_changes, ace_tree_sum$N)

      # store the proportional time spent on each state for this iteration
      mean_state_time[nrow(mean_state_time) + 1,] <- ace_tree_sum[["times"]]["prop",]

      # store the table of state transitions
      total_tr_counts <- append(total_tr_counts, list(ace_tree_sum$Tr))

    }

    # summarise the tables for the most common host links
    table_names <- dir(temp_dir, pattern = "*.tsv", full.names = TRUE)

    most_common_link <- c()

    for (table_file in table_names) {
      link_table <- read.csv(table_file, sep = "\t", header = TRUE)

      most_common_link <-
        c(most_common_link, link_table[c(1), "host_link"])

    }

    most_common_link_df <-
      as.data.frame(table(most_common_link)) %>% arrange(desc(Freq))

    write.table(
      most_common_link_df,
      file = most_common_table,
      row.names = FALSE,
      quote = FALSE,
      sep = '\t'
    )

    # delete the temporary files
    unlink(temp_dir, recursive = TRUE)

    # summarise the ancestral state reconstruction - create a table that holds the mean time spent on each state
    mean_state_time <- colMeans(setNames(mean_state_time, colnames(ace_tree_sum[["times"]])))

    # append the total state changes at the end
    mean_state_time$Avg_state_changes <- mean(total_state_changes)

    write.table(
      mean_state_time,
      file = state_sum_table,
      row.names = FALSE,
      quote = FALSE,
      sep = '\t'
    )

    # summarise the tables of the state transition counts
    mean_state_tr <- apply(laply(total_tr_counts, as.matrix), c(2, 3), mean)

    write.table(
      mean_state_tr,
      file = state_tr_table,
      row.names = FALSE,
      quote = FALSE,
      sep = '\t'
    )
  }

# run the function
transition_analysis(snp_file, treefile, outgroup, pop_meta, all_out_table, most_common_table, state_sum_table,
                    state_tr_table, clock, genome_size, time_interval, cluster)
