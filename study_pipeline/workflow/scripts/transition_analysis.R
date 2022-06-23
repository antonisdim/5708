# Title     : fst
# Objective : Run Nei's (1982) and Weir and Cockerham's Fst
# Created by: Evangelos A. Dimopoulos"
# Created on: 08/04/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
snp_file <- args[1]
treefile <- args[2]
outgroup <- args[3]
pop_meta <- args[4]
out_plot <- args[5]
out_table <- args[6]
mini <- args[7]

# load libs
library(ggplot2)
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
           out_plot,
           out_table,
           mini) {

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

    netowork_phylo <- vcom[-grep("Node", vcom$sample),]

    # add host trait
    sample_meta_df <-
      population_host_metadata("ecoli_population_21_meta.tsv")

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
      netowork_phylo$network_id[match(ecoli_pwsnps_decon$Taxon1, netowork_phylo$sample)]
    ecoli_pwsnps_decon$Taxon2_network <-
      netowork_phylo$network_id[match(ecoli_pwsnps_decon$Taxon2, netowork_phylo$sample)]

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


    if (mini != TRUE) {
      # plot all pairwise distances as histogram
      ecoli_pop21_hist <-
        ggplot(data = ecoli_pwsnps_decon_filtered_clusters,
               aes(x = dist, fill = network_match)) +
        geom_histogram(bins = 200,
                       alpha = 0.8,
                       position = "identity") +
        theme_minimal()  +
        xlab("Pairwise SNP distance") +
        ylab("Frequency") +
        facet_wrap(~ host_link, scales = "free") +
        scale_fill_discrete(name = "Transmission cluster") +
        theme(axis.title = element_text(size = 15)) +
        theme(axis.text.x = element_text(size = 15)) +
        theme(axis.text.y = element_text(size = 15)) +
        theme(strip.text.x = element_text(size = 15, color = "black"))

      ggsave(
        filename = out_plot,
        ecoli_pop21_hist,
        width = 20,
        height = 16,
        units = "in"
      )

    } else {
      ecoli_pwsnps_decon_filtered_clusters_mini <-
        ecoli_pwsnps_decon_filtered_clusters %>%
        filter(dist < 151)

      ecoli_pop21_hist_mini <-
        ggplot(data = ecoli_pwsnps_decon_filtered_clusters_mini,
               aes(x = dist, fill = network_match)) +
        geom_histogram(bins = 200,
                       alpha = 0.8,
                       position = "identity") +
        theme_minimal()  +
        xlab("Pairwise SNP distance") +
        ylab("Frequency") + geom_vline(xintercept = 15) +
        facet_wrap(~ host_link, scales = "free") +
        scale_fill_discrete(name = "Transmission cluster") +
        theme(axis.title = element_text(size = 15)) +
        theme(axis.text.x = element_text(size = 15)) +
        theme(axis.text.y = element_text(size = 15)) +
        theme(strip.text.x = element_text(size = 15, color = "black"))

      ggsave(
        filename = out_plot,
        ecoli_pop21_hist_mini,
        width = 20,
        height = 16,
        units = "in"
      )

    }

    ecoli_pwsnps_decon_filtered_clusters_mini_pairs <-
      ecoli_pwsnps_decon_filtered_clusters %>%
      filter(dist < 16)

    trans_summ <-
      ecoli_pwsnps_decon_filtered_clusters_mini_pairs %>% count(host_link, sort = TRUE)

    write.table(
      trans_summ,
      file = out_table,
      row.names = FALSE,
      quote = FALSE,
      sep = '\t'
    )

  }


tree_obj <- read_tree(treefile, outgroup)

snp_dist_net(snp_file, tree_obj, outgroup, pop_meta, out_plot, out_table, mini)

sample_meta_df <-
  population_host_metadata("ecoli_population_21_meta.tsv")
host_counts <- sample_meta_df %>% count(Trait, sort = TRUE) %>% arrange(n)
subsample_num <- host_counts[2, c('n')]

for (iter in 1:100) {
  out_plot_iter <-
    paste("iter_hist/trans_hist_pop21_mini_iter_", iter, ".pdf", sep = "")
  out_table_iter <-
    paste("iter_tables/trans_table_pop21_iter_", iter, ".tsv", sep = "")
  mini <- TRUE

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
               out_plot_iter,
               out_table_iter,
               mini)

}


table_names <- dir("iter_tables", pattern = "*.tsv", full.names = TRUE)

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
  file = "most_common_link_iter_summary.tsv",
  row.names = FALSE,
  quote = FALSE,
  sep = '\t'
)











#######################################
### Network definition and analyses ###
#######################################

################################################################
### Create network of all samples using SNP threshold of 200 ###
################################################################

################################################
## Create new df containing nodes for network ##
################################################

#allShapes_nodes_200_snps <- sample_meta_df %>%
#  mutate(id = row_number()) %>%
#  select(id, sample, cluster, ST, Country, Trait)

#############################################################
## Extract pairwise SNP comparisons < 16 (edges of network) ##
#############################################################

#ecoli_pwsnps_decon_filtered_clusters_edges_200_snp <-
#  ecoli_pwsnps_decon_filtered_clusters %>%
#  filter(dist < 15) %>%
#  select(Taxon1, Taxon2, dist)

#ecoli_pwsnps_decon_filtered_clusters_edges_200_snp <-
#  ecoli_pwsnps_decon_filtered_clusters_edges_200_snp %>%
#  left_join(allShapes_nodes_200_snps[, c(1, 2)], by = c("Taxon1" = "sample")) %>%
#  dplyr::rename(from.id = id)

#ecoli_pwsnps_decon_filtered_clusters_edges_200_snp <-
#  ecoli_pwsnps_decon_filtered_clusters_edges_200_snp %>%
#  left_join(allShapes_nodes_200_snps[, c(1, 2)], by = c("Taxon2" = "sample")) %>%
#  select(from.id, to.id = id, dist)

###########################
## Create network object ##
###########################

#all_200_snp_routes <-
#  tbl_graph(nodes = allShapes_nodes_200_snps,
#            edges = ecoli_pwsnps_decon_filtered_clusters_edges_200_snp,
#            directed = TRUE)

#all_200_snp_routes_components <- components(all_200_snp_routes)

#allShapes_nodes_200_snps_networked <- allShapes_nodes_200_snps %>%
#  mutate(network_id = all_200_snp_routes_components$membership)

#ecoli_pwsnps_decon_filtered_networked <- ecoli_pwsnps_decon_filtered

#ecoli_pwsnps_decon_filtered_networked$Taxon1_network <-
#  allShapes_nodes_200_snps_networked$network_id[match(
#    ecoli_pwsnps_decon_filtered_networked$Taxon1,
#    allShapes_nodes_200_snps_networked$sample
#  )]

#ecoli_pwsnps_decon_filtered_networked$Taxon2_network <-
#  allShapes_nodes_200_snps_networked$network_id[match(
#    ecoli_pwsnps_decon_filtered_networked$Taxon2,
#    allShapes_nodes_200_snps_networked$sample
#  )]

#plot(ecoli_pwsnps_decon_filtered_networked)
