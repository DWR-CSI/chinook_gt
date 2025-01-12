#!/usr/bin/env Rscript
# Setup --------------
library(tidyverse)
library(data.table)
library(rubias)

# Functions ------------
clean_sample_name <- function(x) {
  x %>%
    str_remove("_$") %>% # Remove trailing underscore
    str_remove("_S\\d+_?$") %>% # Remove _S123_ suffix
    str_remove("_R\\d+_?$") %>% # Remove _R123_ suffix
    str_remove("_L\\d+_?$") # Remove _L123_ suffix
}

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
unks_numgeno <- read_csv(args[1]) %>%
  mutate_if(is.double, as.integer) %>%
  mutate_if(is.logical, as.integer) %>%
  rename_at(vars(ends_with(".1")), ~ str_remove(., "\\.1$")) %>%
  rename_at(vars(ends_with(".2")), ~ str_replace(., "\\.2$", ".1")) %>%
  mutate(sample_type = "mixture", repunit = NA) %>%
  rename(collection = group, indiv = indiv.ID) %>%
  rename_all(~ gsub("-", ".", .)) # dashes not accepted in column names
ref_baseline <- read_csv(args[2]) %>%
  mutate_if(is.double, as.integer) %>%
  rename_all(~ gsub("-", ".", .))
project_name <- args[3]


show_missing_data <- as.logical(args[4])
ots28_info_file <- args[5]
panel_type <- as.character(args[6])
unks_alphageno <- args[7] %>%
  read_tsv() %>%
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.logical, as.character) %>%
  rename_at(vars(ends_with(".1")), ~ str_remove(., "\\.1$")) %>%
  rename_at(vars(ends_with(".2")), ~ str_replace(., "\\.2$", ".1")) %>%
  mutate(sample_type = "mixture", repunit = NA) %>%
  rename(collection = group, indiv = indiv.ID) %>%
  rename_all(~ gsub("-", ".", .))

ots28_missing_threshold <- 0.5 # If less than this much OTS28 data is missing, consider OTS28 data Intermediate instead of uncertain
gsi_missing_threshold <- 0.6 # If more than this much GSI data is missing, consider GSI data invalid
PofZ_threshold <- 0.8 # If the maximum PofZ is less than this, consider the result ambiguous
Spring_PofZ_threshold <- 0.8 # If the maximum PofZ is less than this, consider the result ambiguous
# Parse OTS28 info file ----------------



if (panel_type == "transition") {
  unks <- unks_numgeno
} else if (panel_type == "full") {
  unks <- unks_alphageno
} else {
  stop(paste0("Panel type ", panel_type, " not recognized. Set panel parameter to 'transition' or 'full'."))
}

unks <- unks %>%
  mutate(indiv = clean_sample_name(indiv))

ots28_info <- read_tsv(ots28_info_file) %>%
  mutate(
    indiv = clean_sample_name(indiv),
    RoSA = if_else(RoSA == "Uncertain" & ots28_missing < ots28_missing_threshold, "Intermediate", RoSA)
  ) %>%
  filter(indiv %in% unks$indiv) %>%
  select(indiv, RoSA, ots28_missing)
# Combine unknowns and reference baseline ----------------
unk_match <- unks %>%
  select(any_of(names(ref_baseline))) # Keep only columns (loci) that are in ref_baseline

if (show_missing_data == TRUE) {
  # Add blank/NA data for columns that are in ref_baseline but missing from unk_match
  missing_cols <- setdiff(names(ref_baseline), names(unk_match))
  for (col in missing_cols) {
    unk_match[[col]] <- "ND" # -9 or "ND" is the missing data code
  }
  ref_match <- ref_baseline
  unk_match <- unk_match %>%
    select(names(ref_match))
} else {
  # Remove columns that are in ref_baseline but missing from unk_match
  ref_match <- ref_baseline %>%
    select(any_of(names(unk_match)))
}

unk_match <- unk_match %>%
  mutate(across(everything(), ~ if_else(. == "ND", NA, .))) %>%
  mutate(
    collection = case_when(
      is.na(collection) ~ "ND",
      TRUE ~ collection
    )
  )
for (col in seq(5, ncol(unk_match), 2)) {
  col1 <- col
  col2 <- col + 1
  missing1 <- is.na(unk_match[, col1])
  missing2 <- is.na(unk_match[, col2])
  both_should_be_missing <- missing1 | missing2
  unk_match[both_should_be_missing, col1] <- NA
  unk_match[both_should_be_missing, col2] <- NA
}

chinook_all <- bind_rows(unk_match, ref_match) %>%
  filter(collection != "Coho")


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
if (nrow(unk_match) == 1) {
  all_full_mix_est <- fix_missing_loci(all_full_mix_est)
}


all_full_mix_est_indiv_posteriors <- all_full_mix_est$indiv_posteriors

all_full_mix_results <- all_full_mix_est$indiv_posteriors %>%
  group_by(indiv, mixture_collection, n_miss_loci, n_non_miss_loci, repunit) %>%
  summarize(
    Prob_repunit = sum(PofZ)
  ) %>%
  arrange(indiv, repunit, Prob_repunit)

write_tsv(all_full_mix_est$indiv_posteriors, file = stringr::str_c(project_name, "_full_mix_posteriors.tsv"))


### export all the PofZ for each repunit

repunit_calls <- all_full_mix_results %>%
  group_by(indiv) %>%
  filter(Prob_repunit == max(Prob_repunit)) %>%
  left_join(ots28_info, by = "indiv") %>%
  mutate(
    fraction_missing = n_miss_loci / (n_miss_loci + n_non_miss_loci),
    final_call = case_when(
      fraction_missing > gsi_missing_threshold ~ "Missing Data",
      Prob_repunit < PofZ_threshold ~ "Ambiguous",
      (RoSA == "Early") & (repunit %in% c("fall", "latefall")) ~ "spring",
      (RoSA == "Late") & (repunit == "spring") ~ "Fall",
      (fraction_missing < gsi_missing_threshold) ~ repunit,
      (fraction_missing >= gsi_missing_threshold) & (RoSA == "Late") ~ "Fall / Late Fall",
      (fraction_missing >= gsi_missing_threshold) & (RoSA == "Early") ~ "Spring / Winter",
      TRUE ~ "Assignment Error"
    )
  ) %>%
  select(indiv, RoSA, ots28_missing, n_non_miss_loci, fraction_missing, best_repunit = repunit, prob_repunit = Prob_repunit, final_call)


raw_max_PofZ <- all_full_mix_est_indiv_posteriors %>%
  group_by(indiv, mixture_collection) %>%
  filter(PofZ == max(PofZ))

raw_repunit_probs <- all_full_mix_results %>%
  group_by(indiv) %>%
  filter(Prob_repunit == max(Prob_repunit)) %>%
  left_join(raw_max_PofZ %>% ungroup() %>% select(indiv, collection, PofZ, log_likelihood, z_score), by = "indiv") %>%
  left_join(ots28_info, by = "indiv") %>%
  mutate(
    tributary = case_when(
      (RoSA == "Early") & (repunit %in% c("fall", "latefall")) & (Prob_repunit >= PofZ_threshold) ~ "Feather River Spring",
      (RoSA == "Late") & (repunit == "spring") & (Prob_repunit >= PofZ_threshold) ~ NA_character_,
      repunit == "spring" ~ collection,
      TRUE ~ NA_character_
    ),
    PofZ = case_when(
      (RoSA == "Early") & (repunit %in% c("fall", "latefall")) ~ NA_real_,
      (RoSA == "Late") & (repunit == "spring") ~ NA_real_,
      repunit == "spring" ~ PofZ,
      TRUE ~ NA_real_
    )
  )


mix_results_wide <- all_full_mix_results %>%
  spread(repunit, Prob_repunit) %>%
  group_by(indiv) %>%
  left_join(ots28_info, by = "indiv") %>% # add in the OTS28 info
  left_join(repunit_calls %>% ungroup() %>% select(indiv, fraction_missing, final_call, probability = prob_repunit), by = "indiv") %>% # add in the final calls
  mutate(across(
    matches("spring|winter|fall|late"),
    ~ replace_na(., 0)
  )) %>% # Replace NA with 0
  mutate(across(
    matches("spring|winter|fall|late"),
    ~ if_else(final_call == "Missing Data", NA_real_, .)
  )) %>% # Replace PofZ with NA if final call is "Missing Data"
  mutate(across(
    matches("spring|winter|fall|late"),
    ~ round(., digits = 2)
  )) %>% # Round to 2 decimal places
  mutate(
    GSI_perc_missing = round(fraction_missing * 100, digits = 1),
    ots28_missing = round(ots28_missing, digits = 1),
    probability = if_else((final_call != "Missing Data"), round(probability, digits = 3), NA_real_)
  ) %>%
  select(
    SampleID = indiv,
    RoSA,
    RoSA_perc_missing = ots28_missing,
    GSI_perc_missing,
    Fall = fall,
    Late_fall = latefall,
    Spring = spring,
    Winter = winter,
    final_call,
    probability
  ) %>%
  left_join(raw_repunit_probs %>% select(SampleID = indiv, tributary, trib_PofZ = PofZ), by = "SampleID") %>%
  mutate(
    trib_PofZ = if_else(!(final_call %in% c("Missing Data", "Fall / Late Fall", "Spring / Winter")), round(trib_PofZ, digits = 3), NA_real_)
  )

write_tsv(
  mix_results_wide,
  file = stringr::str_c(project_name, "_summary.tsv")
)
