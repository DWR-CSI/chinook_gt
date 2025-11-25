#!/usr/bin/env Rscript
library("tidyverse")
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
project <- args[2]
rubias_input_file <- args[3]
ckmrsim_input_file <- args[4]

if (mode == 'ckmrsim_rubias') {
  rubias_input <- read_tsv(rubias_input_file)
ckmrsim_input <- read_tsv(ckmrsim_input_file) %>%
  rename(SampleID = D2_indiv, Parent = D1_indiv, PBT_logl_ratio = logl_ratio, PBT_loci = num_loc)
# Select top match for each offspring from filtered CKMRsim results
ckmrsim_input_top = ckmrsim_input %>%
  group_by(SampleID) %>%
  slice_max(PBT_logl_ratio)

combined_rubias_ckmrsim <- rubias_input %>%
  left_join(ckmrsim_input_top, by = "SampleID")
write_tsv(combined_rubias_ckmrsim, file = paste0(project, "_ckmrsim_rubias_summary.tsv"))
}
