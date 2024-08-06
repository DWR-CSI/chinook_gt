#!/usr/bin/env Rscript

# Load required libraries
library(tidyverse)
library(viridis)

# Get command-line arguments passed by Nextflow
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]  # This will be "." (current directory)
output_dir <- args[2]
n_loci <- as.numeric(args[3])

# Function to process a single idxstats file
process_idxstats <- function(file, n_loci) {
  df <- read.table(file, nrows = n_loci, stringsAsFactors = FALSE)
  ind <- str_extract(basename(file), "^[^_]+")
  df$ind <- ind
  return(df)
}

# Main function to process all files and generate plots
main <- function(input_dir, output_dir, n_loci) {
  # List all idxstats files in the input directory
  files <- list.files(input_dir, pattern = "*idxstats", full.names = TRUE)
  
  # Process all files
  idx_list <- map(files, ~process_idxstats(.x, n_loci))
  
  # Combine into a single data frame
  idx_df <- bind_rows(idx_list)
  names(idx_df) <- c("loc", "len", "reads", "unmapd", "ind")
  
  # Summarize and order data
  locs <- idx_df %>% 
    group_by(loc) %>% 
    summarise(total = sum(reads)) %>% 
    arrange(total)
  
  inds <- idx_df %>% 
    group_by(ind) %>% 
    summarise(total = sum(reads)) %>% 
    arrange(total)
  
  idx_df <- idx_df %>%
    mutate(loc = factor(loc, levels = locs$loc),
           ind = factor(ind, levels = inds$ind))
  
  # Generate plots
  p1 <- ggplot(idx_df, aes(loc, ind, fill = reads)) + 
    geom_raster() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    scale_fill_viridis_c(trans = 'log', option = "viridis") +
    ggtitle("Heatmap of reads per locus and individual")
  
  p2 <- ggplot(idx_df, aes(x=ind, y=reads)) +
    geom_col(fill = "lightblue", colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
    xlab("Sample") +
    ylab("Total reads mapped") +
    ggtitle("Reads per individual")
  
  p3 <- ggplot(idx_df, aes(x=loc, y=reads)) +
    geom_col(fill = "goldenrod", colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    xlab("Marker Name, Amplicon Name, or Chromosome") +
    ylab("Total reads mapped") +
    ggtitle("Reads per locus")
  
  # Save plots
  ggsave(file.path(output_dir, "heatmap.pdf"), p1, width = 12, height = 8)
  ggsave(file.path(output_dir, "reads_per_individual.pdf"), p2, width = 12, height = 8)
  ggsave(file.path(output_dir, "reads_per_locus.pdf"), p3, width = 12, height = 8)
  
  # Export data
  matrix_df <- idx_df %>% 
    select(-len, -unmapd) %>% 
    spread(ind, reads)
  write.table(matrix_df, file = file.path(output_dir, "reads_matrix.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
}

# Run the main function
main(input_dir, output_dir, n_loci)