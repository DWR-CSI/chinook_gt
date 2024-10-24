#!/usr/bin/env Rscript
# Setup --------------
library(tidyverse)
library(data.table)
library(rubias)

# Functions ------------
all_na_cols <- function(df) {
  df %>%
    summarise_all(~ all(is.na(.))) %>%
    gather() %>%
    filter(value == TRUE) %>%
    pull(key)
}

# Get command-line arguments passed by Nextflow
args <- commandArgs(trailingOnly = TRUE)
unks <- read_csv(args[1]) %>%
  mutate_if(is.double, as.integer) %>%
  rename_at(vars(ends_with(".1")), ~ str_remove(., "\\.1$")) %>%
  rename_at(vars(ends_with(".2")), ~ str_replace(., "\\.2$", ".1")) %>%
  mutate(sample_type = "mixture", repunit = NA) %>%
  rename(collection = group, indiv = indiv.ID) %>%
  rename_all(~ gsub("-", ".", .)) # dashes not accepted in column names
ref_baseline <- read_csv(args[2]) %>%
  mutate_if(is.double, as.integer) %>%
  rename_all(~ gsub("-", ".", .))
project_name <- args[3]
reporting_groups <- as.integer(args[4])
if (reporting_groups != 2 && reporting_groups != 4) {
  stop("reporting_groups must be 2 or 4")
}

show_missing_data <- as.logical(args[5])
ots28_info_file <- args[6]

ots28_missing_threshold <- 0.5 # If less than this much OTS28 data is missing, consider OTS28 data Intermediate instead of uncertain
gsi_missing_threshold <- 0.6 # If less than this much GSI data is missing, consider GSI data invalid
PofZ_threshold <- 0.7 # If the maximum PofZ is less than this, consider the result ambiguous
# Parse OTS28 info file ----------------
ots28_info <- read_tsv(ots28_info_file) %>%
  mutate(RoSA = if_else(RoSA == "Uncertain" & ots28_missing < ots28_missing_threshold, "Intermediate", RoSA))

# Parse OTS28 info file ----------------
ots28_info <- read_tsv(ots28_info_file) %>%
  mutate(RoSA = if_else(RoSA == "Uncertain" & ots28_missing < ots28_missing_threshold, "Intermediate", RoSA))


# Combine unknowns and reference baseline ----------------
unk_match <- unks %>%
  select(any_of(names(ref_baseline))) # Keep only columns (loci) that are in ref_baseline

if (show_missing_data == TRUE) {
  # Add blank/NA data for columns that are in ref_baseline but missing from unk_match
  missing_cols <- setdiff(names(ref_baseline), names(unk_match))
  for (col in missing_cols) {
    unk_match[[col]] <- as.integer(-9) # -9 is the missing data code
  }
  ref_match <- ref_baseline
  unk_match <- unk_match %>%
    select(names(ref_match))
} else {
  # Remove columns that are in ref_baseline but missing from unk_match
  ref_match <- ref_baseline %>%
    select(any_of(names(unk_match)))
}

chinook_all <- bind_rows(unk_match, ref_match) %>%
  filter(collection != "Coho") %>%
  mutate(across(everything(), ~ if_else(. == -9, NA, .)))
unk_match <- unk_match %>%
  mutate(across(everything(), ~ if_else(. == -9, NA, .)))


# Matching --------------------
matchy_pairs <- close_matching_samples(
  D = chinook_all,
  gen_start_col = 5,
  min_frac_non_miss = 0.85,
  min_frac_matching = 0.94
)

matchy_pairs_out <- matchy_pairs %>%
  arrange(desc(num_non_miss), desc(num_match))

write_tsv(matchy_pairs_out, file = stringr::str_c(project_name, "_matchy_pairs.tsv"))

# Estimate mixtures --------------------
combined_results <- tibble()
if (any(ots28_info$baseline == "SW")) {
  SW_baseline <- ref_match %>%
    filter(
      str_detect(repunit, "Spring") | str_detect(repunit, "Winter")
    )
  SW_unks <- ots28_info %>%
    filter(baseline == "SW") %>%
    pull(indiv)
  SW_unk_match <- unk_match %>%
    filter(indiv %in% SW_unks)
  SW_na_cols <- intersect(all_na_cols(SW_unk_match), all_na_cols(SW_baseline))
  SW_unk_match <- SW_unk_match %>%
    select(-any_of(SW_na_cols))
  SW_baseline <- SW_baseline %>%
    select(-any_of(SW_na_cols))
  SW_mix_est <- infer_mixture(
    reference = SW_baseline,
    mixture = SW_unk_match,
    gen_start_col = 5
  )


  combined_results <- bind_rows(combined_results, SW_mix_est$indiv_posteriors)
}

if (any(ots28_info$baseline == "FLF")) {
  FLF_baseline <- ref_match %>%
    filter(
      str_detect(repunit, "Fall") | str_detect(repunit, "Late")
    )
  FLF_unks <- ots28_info %>%
    filter(baseline == "FLF") %>%
    pull(indiv)
  FLF_unk_match <- unk_match %>%
    filter(indiv %in% FLF_unks)
  # Return column names where all values are NA
  FLF_na_cols <- intersect(all_na_cols(FLF_unk_match), all_na_cols(FLF_baseline))
  FLF_unk_match <- FLF_unk_match %>%
    select(-any_of(FLF_na_cols))
  FLF_baseline <- FLF_baseline %>%
    select(-any_of(FLF_na_cols))

  FLF_mix_est <- infer_mixture(
    reference = FLF_baseline,
    mixture = FLF_unk_match,
    gen_start_col = 5
  )
  combined_results <- bind_rows(combined_results, FLF_mix_est$indiv_posteriors)
}

if (any(ots28_info$baseline == "Full")) {
  full_unks <- ots28_info %>%
    filter(baseline == "Full") %>%
    pull(indiv)
  full_unk_match <- unk_match %>%
    filter(indiv %in% full_unks)
  full_na_cols <- intersect(all_na_cols(full_unk_match), all_na_cols(ref_match))
  full_unk_match <- full_unk_match %>%
    select(-any_of(full_na_cols))
  full_ref_match <- ref_match %>%
    select(-any_of(full_na_cols))
  full_mix_est <- infer_mixture(
    reference = full_ref_match,
    mixture = full_unk_match,
    gen_start_col = 5
  )
  combined_results <- bind_rows(combined_results, full_mix_est$indiv_posteriors)
}

# For troubleshooting only, Full baseline is used here on all samples
all_full_mix_est <- infer_mixture(
  reference = ref_match,
  mixture = unk_match,
  gen_start_col = 5
)

all_full_mix_results <- all_full_mix_est$indiv_posteriors %>%
  group_by(indiv, mixture_collection) %>%
  filter(PofZ == max(PofZ)) %>%
  arrange(indiv, collection, repunit, PofZ)

write_tsv(all_full_mix_results, file = stringr::str_c(project_name, "_full_mix_estimates.tsv"))
