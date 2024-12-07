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
  filter(indiv %in% unks$indiv)
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
      select(-any_of(full_na_cols)) %>%
      mutate(
        repunit = case_when(
          collection == "ColemanLF" ~ "latefall",
          TRUE ~ repunit
        )
      ) %>%
      mutate(collection = repunit)
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
  if (any(ots28_info$baseline == "SW")) {
    SW_baseline <- ref_match %>%
      filter(
        collection %in% c("ButteSp", "FRHsp", "MillDeerSp", "SacWin", "FRHfall")
      ) %>%
      mutate(
        repunit = case_when(
          collection %in% c("FRHsp", "FRHfall") ~ "spring",
          repunit == "winter" ~ "winter",
          TRUE ~ repunit
        ),
        collection = case_when(
          collection == "SacWin" ~ "SacWin",
          TRUE ~ repunit # Set collection to be the same as repunit for other collections
        )
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
        collection %in% c(
          "ButteFall",
          "ColemanLF",
          "FRHfall",
          "FRHsp",
          "MillDeerFall",
          "SanJoaquinFall"
        )
      ) %>%
      mutate(
        repunit = case_when(
          collection == "ColemanLF" ~ "latefall",
          TRUE ~ repunit
        )
      ) %>%
      mutate(
        collection = case_when(
          collection == "ColemanLF" ~ "latefall",
          collection %in% c(
            "ButteFall",
            "ColemanLF",
            "FRHfall",
            "FRHsp",
            "MillDeerFall",
            "SanJoaquinFall"
          ) ~ "fall",
          TRUE ~ repunit
        )
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
        fraction_missing > gsi_missing_threshold ~ "Missing Data",
        PofZ < PofZ_threshold ~ "Ambiguous",
        fraction_missing < gsi_missing_threshold ~ repunit,
        TRUE ~ "Assignment Error"
      )
  )

indivs_too_much_missing_data <- mix_results %>%
  mutate(total_loci = n_miss_loci + n_non_miss_loci) %>% # calculate total number of loci
  filter(n_miss_loci / total_loci > gsi_missing_threshold) %>% # filter for individuals with more than 70% missing data
  pull(indiv) # get the indivs with too much missing data

### export all the PofZ for each repunit

mix_results_long <- combined_results[, c(2:3, 5)]

mix_results_wide <- mix_results_long %>%
  spread(repunit, PofZ) %>%
  group_by(indiv) %>%
  left_join(ots28_info, by = "indiv") %>% # add in the OTS28 info
  left_join(mix_results %>% ungroup() %>% select(indiv, fraction_missing, final_call), by = "indiv") %>% # add in the final calls
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
    ots28_missing = round(ots28_missing, digits = 1)
  ) %>%
  select(
    SampleID = indiv,
    RoSA,
    RoSA_perc_missing = ots28_missing,
    GSI_baseline = baseline,
    GSI_perc_missing,
    Fall = fall,
    Late_fall = latefall,
    Spring = spring,
    Winter = winter,
    final_call
  ) # reorder columns
spring_indivs <- mix_results_wide %>%
  filter(final_call == "spring") %>%
  pull(SampleID)

if ((length(spring_indivs) > 0) && (panel_type == "full")) {
  spring_indiv_data <- unk_match %>%
    filter(indiv %in% spring_indivs)

  # Create Spring trib baseline
  spring_trib_baseline <- ref_match %>%
    filter(
      collection %in% c("ButteSp", "FRHsp", "MillDeerSp")
    ) %>%
    mutate(
      repunit = collection
    )
  spring_trib_na_cols <- intersect(all_na_cols(spring_indiv_data), all_na_cols(spring_trib_baseline))

  spring_trib_baseline <- spring_trib_baseline %>%
    select(-any_of(spring_trib_na_cols))

  spring_indiv_data <- spring_indiv_data %>%
    select(-any_of(spring_trib_na_cols))
  # Run Rubias with Spring unks and Spring trib baseline
  Spring_trib_mix_est <- infer_mixture(
    reference = spring_trib_baseline,
    mixture = spring_indiv_data,
    gen_start_col = 5
  )
  if (nrow(spring_indiv_data) == 1) {
    Spring_trib_mix_est <- fix_missing_loci(Spring_trib_mix_est)
  }
  # Left join the results with mix_results_wide
  Spring_trib_final_calls <- Spring_trib_mix_est$indiv_posteriors %>%
    group_by(indiv, mixture_collection) %>%
    filter(PofZ == max(PofZ)) %>%
    arrange(indiv, collection, repunit, PofZ) %>%
    mutate(
      total_loci = n_miss_loci + n_non_miss_loci,
      fraction_missing = n_miss_loci / total_loci
    ) %>%
    filter(PofZ == max(PofZ)) %>%
    arrange(indiv, collection, repunit, PofZ) %>%
    mutate(
      trib_final_call =
        case_when(
          fraction_missing > gsi_missing_threshold ~ "Missing Data",
          PofZ < Spring_PofZ_threshold ~ "Ambiguous",
          fraction_missing < gsi_missing_threshold ~ repunit,
          TRUE ~ "Assignment Error"
        )
    ) %>%
    ungroup() %>%
    select(SampleID = indiv, trib_final_call)
  Spring_trib_wide <- Spring_trib_mix_est$indiv_posteriors[, c(2:3, 5)] %>%
    rename(SampleID = indiv) %>%
    spread(repunit, PofZ) %>%
    group_by(SampleID) %>%
    left_join(Spring_trib_final_calls, by = "SampleID") %>% # add in the final calls
    mutate(across(matches("ButteSp|FRHsp|MillDeerSp"), ~ round(., digits = 2))) %>% # Replace NA with 0
    select(matches("SampleID|ButteSp|FRHsp|MillDeerSp|trib_final_call"))

  mix_results_wide <- mix_results_wide %>%
    left_join(Spring_trib_wide, by = "SampleID")
}
mix_results_wide <- mix_results_wide %>%
  mutate(
    final_call = case_when(
      GSI_baseline == "FLF" & final_call == "Ambiguous" ~ "Fall / Late Fall",
      TRUE ~ final_call
    )
  )
write_tsv(
  mix_results_wide,
  file = stringr::str_c(project_name, "_summary.tsv")
)
