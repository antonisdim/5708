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

# load libs
library(igraph)
library(ape)
library(dplyr)
library(tidygraph)

# read the utility functions
source("scripts/utilities.R")

snp_dist_net <-
  function(snp_dist,
           tree_obj,
           outgroup,
           pop_meta,
           out_table) {

    # read the snp distance file
    ecoli_pwsnps <- read.csv(snp_dist, sep = '\t', header = TRUE)

    # convert pairsnp output into dataframe
    ecoli_pwsnps_idx <- ecoli_pwsnps
    rownames(ecoli_pwsnps_idx) <- ecoli_pwsnps_idx$X.
    ecoli_pwsnps_idx <-
      ecoli_pwsnps_idx[, -which(names(ecoli_pwsnps_idx) %in% c("X", "X."))]

    ecoli_pwsnps_decon <-
      data.frame(t(combn(names(ecoli_pwsnps_idx), 2)), dist = t(ecoli_pwsnps_idx)[lower.tri(ecoli_pwsnps_idx)])
    colnames(ecoli_pwsnps_decon) <- c("Taxon1", "Taxon2", "dist")
    ecoli_pwsnps_decon$Taxon1 <-
      gsub('# ', "", ecoli_pwsnps_decon$Taxon1)

    # read the corresponding tree, create a network out of it and then convert it to df
    tree_no_out <- tree_obj

    tree_no_out$node.label <- NULL
    xg <- as.igraph(tree_no_out, directed = FALSE)
    comm <- cluster_fast_greedy(xg)
    V(xg)$community <- membership(comm)

    vcom <- data.frame(as_ids(V(xg)), V(xg)$community)
    colnames(vcom) <- c("sample", "network_id")

    network_phylo <- vcom[-grep("Node", vcom$sample),]

    # add host trait
    sample_meta_df <-
      population_host_metadata(pop_meta)

    # match sample and host
    ecoli_pwsnps_decon$Taxon1_host <-
      sample_meta_df$Trait[match(ecoli_pwsnps_decon$Taxon1, sample_meta_df$sample)]
    ecoli_pwsnps_decon$Taxon2_host <-
      sample_meta_df$Trait[match(ecoli_pwsnps_decon$Taxon2, sample_meta_df$sample)]

    # match sample and trial area
    ecoli_pwsnps_decon$Taxon1_area <-
      sample_meta_df$Country[match(ecoli_pwsnps_decon$Taxon1, sample_meta_df$sample)]
    ecoli_pwsnps_decon$Taxon2_area <-
      sample_meta_df$Country[match(ecoli_pwsnps_decon$Taxon2, sample_meta_df$sample)]

    # match sample and network_id
    ecoli_pwsnps_decon$Taxon1_network <-
      network_phylo$network_id[match(ecoli_pwsnps_decon$Taxon1, network_phylo$sample)]
    ecoli_pwsnps_decon$Taxon2_network <-
      network_phylo$network_id[match(ecoli_pwsnps_decon$Taxon2, network_phylo$sample)]

    # create new df with columns for area_match and host_match to help with separating and filtering data for plots
    ecoli_pwsnps_decon_filtered <- ecoli_pwsnps_decon %>%
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

    # create new df filtering out samples that aren't in the same trial area                                                  ##
    ecoli_pwsnps_decon_filtered_clusters <-
      ecoli_pwsnps_decon_filtered %>%
      filter(!is.na(Taxon1_host)) %>%
      filter(!is.na(Taxon2_host)) %>%
      filter(!is.na(network_match)) %>%
      filter(!grepl("Other", host_link))

    ecoli_pwsnps_decon_filtered_clusters_mini_pairs <-
      ecoli_pwsnps_decon_filtered_clusters %>%
      filter(dist < 16)

    trans_summary <-
      ecoli_pwsnps_decon_filtered_clusters_mini_pairs %>% count(host_link, sort = TRUE)

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
           most_common_table) {

    tree_obj <- read_tree(treefile, outgroup)
    snp_dist_net(snp_file, tree_obj, outgroup, pop_meta, all_out_table)

    sample_meta_df <- population_host_metadata(pop_meta)
    host_counts <- sample_meta_df %>% count(Trait, sort = TRUE) %>% arrange(n)
    subsample_num <- host_counts[2, c('n')]

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

      for (host in unique(sample_meta_df$Trait)) {
        host_count <- nrow(sample_meta_df[sample_meta_df$Trait == host, ])

        if (host_count < subsample_num) {
          nsample = host_count
        } else {
          nsample = subsample_num
        }

        host_subsample <-
          sample_meta_df %>% filter(Trait == host) %>% sample_n(., nsample)
        subsample_list[[host]] <- host_subsample

      }

      subsample_df <- do.call(rbind, subsample_list)

      to_drop <-
        tree_obj$tip.label[!tree_obj$tip.label %in% subsample_df$sample]
      subset_tree <- tree_obj

      for (tip in to_drop) {
        subset_tree <-
          drop.tip(subset_tree,
                   c(tip),
                   rooted = TRUE,
                   collapse.singles = TRUE)


      }

      snp_dist_net(snp_file,
                   subset_tree,
                   outgroup,
                   pop_meta,
                   out_table_iter)

    }

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

    unlink(temp_dir, recursive=TRUE)

  }

# run the install function
transition_analysis(snp_file, treefile, outgroup, pop_meta, all_out_table, most_common_table)
