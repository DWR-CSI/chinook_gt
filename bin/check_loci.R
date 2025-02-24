#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)
library(viridis)
library(cowplot)

# Constants
TOTAL_PANEL_LOCI <- 204
MIN_READS <- 10
POOR_PERFORMANCE_THRESHOLD <- 0.5

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

# Read and process data
data <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE) %>%
  filter(loc != '*')

# Verify locus count
if (nrow(data) != TOTAL_PANEL_LOCI) {
  warning(sprintf("Expected %d loci but found %d in input file", 
                 TOTAL_PANEL_LOCI, nrow(data)))
}

# Convert to long format
long_data <- data %>%
  pivot_longer(cols = -loc, 
              names_to = "sample", 
              values_to = "reads")

# Calculate sample-level metrics
sample_stats <- long_data %>%
  group_by(sample) %>%
  summarize(
    total_reads = sum(reads),
    mean_reads = mean(reads),
    median_reads = median(reads),
    sd_reads = sd(reads),
    cv = sd_reads / mean_reads * 100,
    non_zero_loci = sum(reads > MIN_READS),
    total_loci = TOTAL_PANEL_LOCI,
    pct_success = (non_zero_loci / TOTAL_PANEL_LOCI) * 100,
    shannon = diversity(reads, index = "shannon"),
    pielou = shannon / log(TOTAL_PANEL_LOCI),
    simpson = diversity(reads, index = "simpson")
  ) %>%
  arrange(desc(total_reads))

# Calculate locus-level metrics
locus_stats <- long_data %>%
  group_by(loc) %>%
  summarize(
    total_reads = sum(reads),
    mean_reads = mean(reads),
    median_reads = median(reads),
    sd_reads = sd(reads),
    cv = sd_reads / mean_reads * 100,
    non_zero_samples = sum(reads > MIN_READS),
    total_samples = n(),
    pct_success = (non_zero_samples / total_samples) * 100
  ) %>%
  arrange(desc(pct_success))

# Identify problematic loci
problem_loci <- locus_stats %>%
  filter(pct_success < (POOR_PERFORMANCE_THRESHOLD * 100)) %>%
  arrange(pct_success)

# Create visualizations
p1 <- ggplot(sample_stats, aes(x = reorder(sample, -total_reads), y = total_reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total Reads per Sample",
       x = "Sample",
       y = "Total Reads")

p2 <- ggplot(sample_stats, aes(x = total_reads, y = pct_success)) +
  geom_point(aes(size = non_zero_loci), alpha = 0.6) +
  geom_smooth(method = "lm", alpha = 0.2) +
  theme_minimal() +
  labs(title = "Success Rate vs Total Reads",
       x = "Total Reads",
       y = "Percent Success",
       size = "Successful Loci")

# Create read depth heatmap
heatmap_data <- long_data %>%
  mutate(read_category = cut(reads, 
                            breaks = c(-Inf, 0, 10, 100, 1000, Inf),
                            labels = c("0", "1-10", "11-100", "101-1000", ">1000")))

p3 <- ggplot(heatmap_data, 
             aes(x = reorder(sample, -reads), 
                 y = reorder(loc, reads), 
                 fill = read_category)) +
  geom_tile() +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8)) +
  labs(title = "Read Depth Heatmap",
       x = "Sample",
       y = "Locus",
       fill = "Read Count")

p4 <- ggplot(sample_stats, aes(x = reorder(sample, -pielou), y = pielou)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Amplification Evenness by Sample",
       x = "Sample",
       y = "Pielou's Evenness")

# Save results
write.csv(sample_stats, 
          file.path(output_dir, "sample_statistics.csv"), 
          row.names = FALSE)
write.csv(locus_stats, 
          file.path(output_dir, "locus_statistics.csv"), 
          row.names = FALSE)
write.csv(problem_loci, 
          file.path(output_dir, "problematic_loci.csv"), 
          row.names = FALSE)

# Save plots
plots <- plot_grid(p1, p2, p3, p4, ncol = 2)
ggsave(file.path(output_dir, "loci_qc_plots.pdf"), 
       plots, 
       width = 15, 
       height = 12)

# Generate summary report
sink(file.path(output_dir, "loci_qc_summary.txt"))
cat("Amplicon Sequencing QC Analysis Summary\n")
cat("======================================\n\n")
cat(sprintf("Total Samples: %d\n", nrow(sample_stats)))
cat(sprintf("Total Loci: %d\n", TOTAL_PANEL_LOCI))
cat(sprintf("Mean Reads per Sample: %.0f\n", mean(sample_stats$total_reads)))
cat(sprintf("Median Reads per Sample: %.0f\n", median(sample_stats$total_reads)))
cat(sprintf("Mean Success Rate: %.1f%%\n", mean(sample_stats$pct_success)))
cat(sprintf("Number of Problematic Loci: %d\n", nrow(problem_loci)))
cat(sprintf("Mean Amplification Evenness: %.3f\n", mean(sample_stats$pielou)))
sink()