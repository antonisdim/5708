# Title     : gene_carriage_comparison
# Objective : Do a first pass on the N of amr or vf genes that are carried by isolates each year
# Created by: Evangelos A. Dimopoulos"
# Created on: 16/11/2023

# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
pop_metadata <- args[1]
gene_table <- args[2]
dataset <- args[3]

# load libs
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)

# read the utility functions
source("/Users/ed601/5708/study_pipeline/workflow/scripts/utilities.R")

theme_Publication <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.key.width = unit(1.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

c36 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", 
  "#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089", 
  "#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#76716E", 
  "#DEB478", "#976533", "#D3D5D0","#DD3C51", "#B40F20"
)


# define a basename for the outputs
outfile_base <- tools::file_path_sans_ext(basename(gene_table))

# start producing produce the area curves for the dataset

gene_df <- read.csv(gene_table, sep='\t', header=TRUE)
gene_counts <- gene_df[, c('X.FILE' ,'NUM_FOUND')]

pop_meta <- population_host_metadata(pop_metadata)
pop_meta <- pop_meta %>% right_join(gene_counts, by=join_by('sample' == 'X.FILE'))
if (dataset == 'sb27') pop_meta <- pop_meta[pop_meta$Dataset == 'SB27',]


# group by per year and baps cluster and plot

data <- pop_meta  %>%
  group_by(Collection_Year, cluster) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n / sum(n))*100, cluster = as.character(cluster))

data$cluster <- factor(data$cluster , levels=c("1", "2", "3", "4", "5", "6",
                                               "7", "8", "9", "10", "11", "12",
                                               "13", "14", "15", "16", "17", 
                                               "18", "19", "20", "21", "22", 
                                               "23", "24", "25", "26", "27",
                                               "28", "29", "30", "31", "32",
                                               "33", "34", "35", "36"))

# produce and save the plot

area_curve_plot <- ggplot(data[(data$Collection_Year>=2017) & 
              (data$Collection_Year<=2021),], 
       aes(x=Collection_Year, y=percentage, fill=cluster)) + 
  geom_area(alpha=0.8 , size=1, colour="black") +
  scale_fill_manual(values=c36) +
  scale_x_continuous(breaks=seq(2016, 2022, 1)) +
  theme_Publication()

ggsave(paste(outfile_base, '_area_curves.pdf', sep=''), plot=area_curve_plot, width=9, height=9)


# calculate the trends of gene content per year 

trend <- spread(subset(data, select=-c(n)), Collection_Year, percentage)
#trend <- spread(subset(data, select=-c(n, amr_count)), Collection_Year, amr_percentage)

max_st <- pop_meta %>% 
  group_by(cluster) %>%
  count(ST) %>%
  slice(which.max(n))

max_host <- pop_meta %>% 
  group_by(cluster) %>%
  count(Trait) %>%
  slice(which.max(n))

annot_df <- max_st %>% right_join(max_host, by=('cluster')) %>% 
  mutate(cluster = as.character(cluster)) %>%
  right_join(trend, by=('cluster'))

write.table(annot_df, file = paste(outfile_base, '_trend.tsv', sep=''), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


# AMR or VF gene barplots - in a for loop 

for (source in c("All samples", "Human", "Animal")) {
  
  pop_freq <- pop_meta
  
  if (source != 'All samples') pop_freq <- pop_meta[pop_meta$Source == source,]
  
  # group by per year and baps cluster and plot
  gene_freq_data <- pop_freq  %>%
    group_by(Collection_Year) %>%
    summarise(n = n(), gene_count = sum(NUM_FOUND)) %>%
    mutate(gene_percentage = (gene_count / sum(gene_count))*100)
      
      
    bar_plot <- ggplot(data=gene_freq_data[(gene_freq_data$Collection_Year>=2016) & 
                                             (gene_freq_data$Collection_Year<=2021),], 
                       aes(x=Collection_Year, y=gene_percentage)) + 
      geom_bar(stat='identity', fill="skyblue2") + 
      ggtitle(source) + 
      ylab('Gene carriage percentage') + 
      xlab("Collection year") + 
      theme_Publication()
  
    if (source == 'All samples') out_fix <- 'all_samples' else out_fix <- source
  
  ggsave(paste(outfile_base, '_', source, '_freq_barplot.pdf', sep=''), plot=bar_plot, width=8, height=8)
  
}


# AMR or VF gene barplots - for each baps cluster separately 

freqs_per_cluster <- pop_meta  %>%
  group_by(Collection_Year, cluster) %>%
  summarise(n = n(), gene_count = sum(NUM_FOUND)) %>%
  mutate(test = sum(gene_count), gene_percentage = (gene_count / sum(gene_count))*100)

freqs_per_cluster_plot <- ggplot(data=freqs_per_cluster[(freqs_per_cluster$Collection_Year>=2016) & 
                       (freqs_per_cluster$Collection_Year<=2021),], 
       aes(x=Collection_Year, y=gene_percentage)) + 
  geom_bar(stat='identity', fill="skyblue2") +
  facet_wrap(~cluster) +
  ggtitle("Gene carriage per baps cluster") +
  ylab('Gene carriage percentage') +
  xlab("Collection year") +
  theme_Publication() 

ggsave(paste(outfile_base, '_freq_per_baps_barplot.pdf', sep=''), plot=freqs_per_cluster_plot, width=12, height=10)


# Use a normalized estimate of the gene carriage 

# substitute '.' with NA 
gene_df[gene_df == "."] <- NA

# melt the dataframe 
gene_df_melt <- melt(gene_df, na.rm = FALSE, value.name = "gene_id_percent", 
                     variable.name = "gene", id = c("X.FILE", "NUM_FOUND"))
gene_df_melt <- gene_df_melt %>% 
  right_join(pop_meta[, c("sample", "Collection_Year")], by=join_by('X.FILE' == 'sample'))

# get the total number of unique genes in the dataset 
unique_genes_total <- length(unique(gene_df_melt$gene))

# define the dataframe tos tore the ratios
gene_ratio_df <- data.frame(year = double(), gene_carriage_ratio = double())
                         
# iterate through it per year 
for (year in sort(unique(gene_df_melt$Collection_Year))) {
  per_year <- gene_df_melt[gene_df_melt$Collection_Year == year,]
  
  # exclude rows/individuals where the gene does not exist - otherwise unique_genes_yearly == unique_genes_total always 
  per_year <- per_year %>% drop_na(gene_id_percent)
  
  # count the number of unique genes in a year 
  unique_genes_yearly <- length(unique(per_year$gene))
  
  # count the number of individuals in that year
  unique_indiv_year <- length(unique(per_year$X.FILE))
  
  # calculate ratio number of unique genes in a year / number of unique genes in dataset / number of isolates in a year
  gene_ratio <- (unique_genes_yearly/unique_genes_total)/unique_indiv_year
  
  print(gene_ratio)
  
  # store the results 
  res_vector <- c(year, gene_ratio)
  gene_ratio_df[nrow(gene_ratio_df) + 1,] <- res_vector
}

# append to the p-values dataframe
write.table(gene_ratio_df, file = paste(outfile_base, '_ratios.tsv', sep=''), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


