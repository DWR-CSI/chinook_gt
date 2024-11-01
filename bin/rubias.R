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

fix_missing_loci <- function(mix_est) {
  # Get unique loci values and names
  if (is.matrix(mix_est$indiv_posteriors$missing_loci)) {
    # For each repunit/collection combination, create a single row with list of missing loci
    fixed_posteriors <- mix_est$indiv_posteriors %>%
      group_by(repunit, collection) %>%
      slice(1) %>% # Take just the first row for each repunit/collection combo
      ungroup()

    # Convert the missing_loci column to a list
    fixed_posteriors$missing_loci <- list(
      structure(mix_est$indiv_posteriors$missing_loci[, 1],
        names = rownames(mix_est$indiv_posteriors$missing_loci)
      )
    )

    # Replace the indiv_posteriors in the original object
    mix_est$indiv_posteriors <- fixed_posteriors
  }
  return(mix_est)
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
panel_type <- as.character(args[7])

ots28_missing_threshold <- 0.5 # If less than this much OTS28 data is missing, consider OTS28 data Intermediate instead of uncertain
gsi_missing_threshold <- 0.6 # If more than this much GSI data is missing, consider GSI data invalid
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
combined_results <- tibble(
  mixture_collection = character(),
  indiv = character(),
  repunit = character(),
  collection = character(),
  PofZ = double(),
  log_likelihood = double(),
  z_score = double(),
  n_non_miss_loci = integer(),
  n_miss_loci = integer(),
  missing_loci = list()
)

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

if (panel_type == "transition") {
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

    if (nrow(SW_unk_match) == 1) {
      SW_mix_est <- fix_missing_loci(SW_mix_est)
    }
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
    if (nrow(FLF_unk_match) == 1) {
      FLF_mix_est <- fix_missing_loci(FLF_mix_est)
    }
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
    if (nrow(full_unk_match) == 1) {
      full_mix_est <- fix_missing_loci(full_mix_est)
    }
    combined_results <- bind_rows(combined_results, full_mix_est$indiv_posteriors)
  }
} else if (panel_type == "full") {
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
    if (nrow(full_unk_match) == 1) {
      full_mix_est <- fix_missing_loci(full_mix_est)
    }
    combined_results <- bind_rows(combined_results, full_mix_est$indiv_posteriors)
  }
} else {
  stop(paste0("Panel type ", panel_type, " not recognized. Set panel parameter to 'transition' or 'full'."))
}

mix_results <- combined_results %>%
  group_by(indiv, mixture_collection) %>%
  mutate(
    total_loci = n_miss_loci + n_non_miss_loci,
    fraction_missing = n_miss_loci / total_loci
  ) %>%
  filter(PofZ == max(PofZ)) %>%
  arrange(indiv, collection, repunit, PofZ) %>%
  mutate(
    final_call =
      case_when(
        PofZ < PofZ_threshold ~ "Ambiguous",
        fraction_missing < gsi_missing_threshold ~ repunit,
        TRUE ~ "Missing Data"
      )
  )

indivs_too_much_missing_data <- mix_results %>%
  mutate(total_loci = n_miss_loci + n_non_miss_loci) %>% # calculate total number of loci
  filter(n_miss_loci / total_loci > gsi_missing_threshold) %>% # filter for individuals with more than 70% missing data
  pull(indiv) # get the indivs with too much missing data




### export all the PofZ for each repunit

mix_results_long <- combined_results[, c(2:3, 5)]
final_calls <- mix_results %>%
  ungroup() %>%
  select(indiv, repunit, n_non_miss_loci, n_miss_loci, total_loci, fraction_missing, final_call)

mix_results_wide <- mix_results_long %>%
  spread(repunit, PofZ) %>%
  left_join(ots28_info, by = "indiv") %>% # add in the OTS28 info
  left_join(mix_results %>% ungroup() %>% select(indiv, fraction_missing, final_call), by = "indiv") %>% # add in the final calls
  mutate(across(
    matches("Spring|Winter|Fall|Late"),
    ~ replace_na(., 0)
  )) %>% # Replace NA with 0
  mutate(across(
    matches("Spring|Winter|Fall|Late"),
    ~ if_else(final_call == "Missing Data", NA_real_, .)
  )) %>% # Replace PofZ with NA if final call is "Missing Data"
  mutate(across(
    matches("Spring|Winter|Fall|Late"),
    ~ round(., digits = 2)
  )) %>% # Round to 2 decimal places
  mutate(
    GSI_perc_missing = fraction_missing * 100
  ) %>%
  select(
    SampleID = indiv,
    RoSA,
    RoSA_perc_missing = ots28_missing,
    GSI_baseline = baseline,
    GSI_perc_missing,
    Fall = CV_Fall,
    Late_fall = CV_Late_fall,
    Spring = CV_Spring,
    Winter = CV_Winter,
    final_call
  ) # reorder columns
write_tsv(
  mix_results_wide,
  file = stringr::str_c(project_name, "_summary.tsv")
)
