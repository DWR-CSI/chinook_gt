#!/usr/bin/env Rscript
# Setup --------------
library(tidyverse)

# Functions ------------

# Function to extract summary statistics
extract_summary_stats <- function(lines) {
  # Find the line with "Proportion of membership of each pre-defined"
  start_line <- which(str_detect(lines, "Proportion of membership of each pre-defined"))[1] + 6
  end_line <- start_line + 4

  # Extract the relevant lines
  summary_lines <- lines[start_line:end_line]

  # Parse the summary statistics
  summary_stats <- read.table(
    text = summary_lines, header = FALSE,
    col.names = c("Pop", "Cluster1", "Cluster2", "NumIndividuals")
  ) %>%
    mutate(
      Pop = as.numeric(str_remove(Pop, ":")),
      NumIndividuals = as.numeric(NumIndividuals)
    ) %>%
    mutate(Cluster_diff = (Cluster1 - Cluster2))
  # groups 1 and 2 are late, 3 and 4 are early
  return(summary_stats)
}

# Function to extract inferred ancestry data
extract_ancestry_data <- function(lines) {
  # Find the line with "Inferred ancestry of individuals:"
  start_line <- which(str_detect(lines, "Inferred ancestry of individuals:"))[1] + 2
  end_line <- which(str_detect(lines, "Estimated Allele Frequencies in each cluster"))[1] - 3

  # Extract the relevant lines
  ancestry_lines <- lines[start_line:end_line]

  # Parse the ancestry data
  ancestry_data <- read.table(text = ancestry_lines, header = FALSE, col.names = c("Index", "Label", "Missing", "Pop", "sep", "Cluster1", "Cluster2")) %>%
    mutate(
      Missing = as.numeric(str_remove_all(Missing, "[()]")),
      Pop = as.numeric(str_remove(Pop, ":"))
    ) %>%
    select(-sep)

  return(ancestry_data)
}

which_RoSA_informed_baseline <- function(RoSA) {
  if (RoSA == "Late") {
    return("FLF")
    # needs to continue on to use baseline that only includes Fall + Late Fall pops
  } else if (RoSA == "Early") {
    return("SW")
  } else if (RoSA == "Intermediate") {
    return("Full")
    # Does not need to continue onto another baseline
  } else if (RoSA == "Uncertain") {
    return("Full")
  } else {
    return("Error")
  }
}

# Main ------------
# Get command-line arguments passed by Nextflow
args <- commandArgs(trailingOnly = TRUE)

structure_output_path <- args[1]
project_name <- args[2]
structure_input_path <- args[3]
missing_threshold <- args[4] * 100 # Percentage. Samples with more missing marked as uncertain.
ots28_output_filename <- stringr::str_c(project_name, "_ots28_report.tsv")
lines <- readLines(structure_output_path)



# Use the function

structure_output <- readLines(structure_output_path)

summary_data <- extract_summary_stats(structure_output)
ancestry_data <- extract_ancestry_data(structure_output)

refpop1_diff <- summary_data %>%
  filter(Pop == 1) %>%
  select(Cluster_diff) %>%
  pull()

refpop2_diff <- summary_data %>%
  filter(Pop == 2) %>%
  select(Cluster_diff) %>%
  pull()

refpop3_diff <- summary_data %>%
  filter(Pop == 3) %>%
  select(Cluster_diff) %>%
  pull()

refpop4_diff <- summary_data %>%
  filter(Pop == 4) %>%
  select(Cluster_diff) %>%
  pull()

if (
  refpop1_diff < -0.7 && refpop2_diff < -0.7 && refpop3_diff > 0.7 && refpop4_diff > 0.7
) {
  ancestry_data <- ancestry_data %>%
    mutate(Ancestry = case_when( # Cluster 1 is early, cluster 2 is late
      Cluster1 > Cluster2 * 5 ~ "Early", # equivalent to Cluster 1 > 83% and Cluster 2 < 17%
      Cluster1 * 5 < Cluster2 ~ "Late",
      Missing < missing_threshold ~ "Intermediate",
      TRUE ~ "Uncertain"
    ))
} else if (refpop1_diff > 0.7 && refpop2_diff > 0.7 && refpop3_diff < -0.7 && refpop4_diff < -0.7) {
  ancestry_data <- ancestry_data %>%
    mutate(Ancestry = case_when( # Cluster 1 is late, cluster 2 is early
      Cluster1 * 5 < Cluster2 ~ "Early",
      Cluster1 > Cluster2 * 5 ~ "Late",
      Missing < missing_threshold ~ "Intermediate",
      TRUE ~ "Uncertain"
    ))
} else {
  ancestry_data <- ancestry_data %>%
    mutate(Ancestry = "Uncertain")
  warning("Unable to determine clustering of the reference populations. All samples will be marked as uncertain.")
}

ancestry_data <- ancestry_data %>%
  mutate(
    Ancestry = ifelse(Missing > missing_threshold, "Uncertain", Ancestry)
  ) %>%
  rename(
    "indiv" = "Label",
    "ots28_missing" = "Missing",
    RoSA = "Ancestry"
  ) %>%
  select(indiv, RoSA, ots28_missing) %>%
  mutate(baseline = map_chr(RoSA, which_RoSA_informed_baseline))


# Replace labels with input labels from original data because Structure truncates the labels.

original_labels <- read_tsv(file = structure_input_path, skip = 1, col_names = FALSE) %>%
  pull(1)
if (length(original_labels) == length(ancestry_data$indiv)) {
  ancestry_data$indiv <- original_labels
} else {
  warning("The number of samples in the ancestry data does not match the number of samples in the original data.")
}

write_tsv(ancestry_data, ots28_output_filename)
