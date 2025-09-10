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

# Chinook species checker function
# Determines species classification based on diagnostic locus genotype
# This locus is diagnostic for Chinook salmon species identification
# 
# Parameters:
#   data: Data frame containing genetic data with diagnostic locus columns
#   diagnostic_locus: Base name of the diagnostic locus (default: "OkiOts_120255.113")
#
# Returns:
#   Data frame with SampleID and Species columns
#   Species classifications:
#     - "Chinook": Both alleles are "GA" (homozygous for Chinook-specific allele)
#     - "Unconfirmed": Either allele is missing data ("ND")  
#     - "non-Chinook": Both alleles are present but at least one is not "GA"
#     - "Error": Catch-all for unexpected cases
chinook_species_checker <- function(data, diagnostic_locus = "OkiOts_120255.113", suffix = ".1") {
  # Check if required columns exist
  allele1_col <- diagnostic_locus
  allele2_col <- paste0(diagnostic_locus, suffix)
  
  if (!all(c(allele1_col, allele2_col) %in% names(data))) {
    warning(paste("Species diagnostic columns", allele1_col, "and/or", allele2_col, "not found. Returning all samples as 'Unconfirmed'"))
    return(data %>% 
           select(SampleID = indiv) %>%
           mutate(Species = "Unconfirmed"))
  }
  
  # Check if indiv column exists
  if (!"indiv" %in% names(data)) {
    stop("Column 'indiv' not found in data")
  }
  
  data %>%
    mutate(
      Species = case_when(
        .data[[allele1_col]] == "GA" & .data[[allele2_col]] == "GA" ~ "Chinook",
        .data[[allele1_col]] == "ND" | .data[[allele2_col]] == "ND" ~ "Unconfirmed",
        is.na(.data[[allele1_col]]) | is.na(.data[[allele2_col]]) ~ "Unconfirmed",
        .data[[allele1_col]] != "GA" | .data[[allele2_col]] != "GA" ~ "non-Chinook",
        TRUE ~ "Error"
      )
    ) %>%
    select(SampleID = indiv, Species)
}

LFAR_checker <- function(data, diagnostic_loci = c("NC_037130.1:864908.865208", "NC_037130.1:1062935.1063235"), suffix = ".1") {
  allele1_cols <- diagnostic_loci
  allele2_cols <- paste0(diagnostic_loci, suffix)
  all_cols <- c(allele1_cols, allele2_cols)
  
  # Ensure all columns exist
  missing_cols <- setdiff(all_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing LFAR check columns: ", paste(missing_cols, collapse = ", "))
  }
  
  data %>%
    mutate(
      LFAR_markers_present = if_else(
        rowSums(!is.na(across(all_of(all_cols))) & across(all_of(all_cols)) != "ND") > 0,
        TRUE, FALSE
      )
    ) %>%
    select(SampleID = indiv, LFAR_markers_present)
}

calculate_heterozygosity <- function(data, gen_start_col = 5) {
  # Calculate heterozygosity for each individual
  het_data <- data %>%
    select(indiv, everything()) %>%
    mutate(
      heterozygosity = map_dbl(1:n(), function(row_idx) {
        # Get genetic data for this individual starting from gen_start_col
        genetic_data <- as.character(unlist(.[row_idx, gen_start_col:ncol(.)]))
        
        # Group into pairs (every 2 columns is a locus)
        n_loci <- length(genetic_data) %/% 2
        valid_loci <- 0
        het_loci <- 0
        
        for (locus in 1:n_loci) {
          allele1_idx <- (locus - 1) * 2 + 1
          allele2_idx <- (locus - 1) * 2 + 2
          
          if (allele1_idx <= length(genetic_data) && allele2_idx <= length(genetic_data)) {
            allele1 <- genetic_data[allele1_idx]
            allele2 <- genetic_data[allele2_idx]
            
            # Check if both alleles are present (not NA or ND)
            if (!is.na(allele1) && !is.na(allele2) && allele1 != "ND" && allele2 != "ND") {
              valid_loci <- valid_loci + 1
              if (allele1 != allele2) {
                het_loci <- het_loci + 1
              }
            }
          }
        }
        
        # Return proportion of heterozygous loci
        if (valid_loci > 0) {
          return(het_loci / valid_loci)
        } else {
          return(NA_real_)
        }
      })
    ) %>%
    select(SampleID = indiv, heterozygosity) %>%
    mutate(heterozygosity = round(heterozygosity, digits = 3))
  
  return(het_data)
}


# Get command-line arguments passed by Nextflow
args <- commandArgs(trailingOnly = TRUE)
ref_baseline <- read_csv(args[1]) %>%
  mutate_if(is.double, as.integer) %>%
  rename_all(~ gsub("-", ".", .))
project_name <- args[2]


show_missing_data <- as.logical(args[3])
ots28_info_file <- args[4]
panel_type <- as.character(args[5])
unks_alphageno <- args[6] %>%
  read_csv() %>%
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.logical, as.character) %>%
  mutate(across(everything(), ~ if_else(. == "NA", "ND", .))) %>%
  rename_at(vars(ends_with(".1")), ~ str_remove(., "\\.1$")) %>%
  rename_at(vars(ends_with(".2")), ~ str_replace(., "\\.2$", ".1")) %>%
  mutate(sample_type = "mixture", repunit = NA) %>%
  rename(collection = group, indiv = indiv.ID) %>%
  rename_all(~ gsub("-", ".", .))

ots28_missing_threshold <- as.numeric(args[7]) * 100 # If less than this much OTS28 data is missing, consider OTS28 data Intermediate instead of uncertain. Multiplied by 100 to get percentage
gsi_missing_threshold <- as.numeric(args[8]) # If more than this much GSI data is missing, consider GSI data invalid
PofZ_threshold <- as.numeric(args[9]) # If the maximum PofZ is less than this, consider the result ambiguous
Spring_PofZ_threshold <- PofZ_threshold # Not currently used, but could be used to set a minimum PofZ for spring trib calls
# Parse OTS28 info file ----------------



if (panel_type == "full") {
  unks <- unks_alphageno
} else {
  stop(paste0("Panel type ", panel_type, " not recognized. Set panel parameter to 'full'."))
}

unks <- unks %>%
  mutate(indiv = clean_sample_name(indiv))


# Species identification and genetic quality assessment ----------------
# Run species checker on unknown samples after sample name cleanup
# This identifies potential non-Chinook samples based on diagnostic locus OkiOts_120255-113
species_results <- chinook_species_checker(unks)

# Calculate genome-wide heterozygosity for all unknown samples
# This provides a measure of genetic diversity and can help identify potential issues including non-Chinook samples
heterozygosity_results <- calculate_heterozygosity(unks)

LFAR_results <- LFAR_checker(unks)

ots28_info <- read_tsv(ots28_info_file) %>%
  mutate(
    indiv = clean_sample_name(indiv),
    RoSA = case_when(
      ots28_missing >= ots28_missing_threshold ~ "Uncertain",
      RoSA == "Uncertain" & ots28_missing < ots28_missing_threshold ~ "Intermediate",
      TRUE ~ RoSA
    )
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
      (fraction_missing > gsi_missing_threshold) & (RoSA == "Late") ~ "Fall / Late Fall",
      (fraction_missing > gsi_missing_threshold) & (RoSA == "Early") ~ "Spring / Winter",
      fraction_missing > gsi_missing_threshold ~ "Missing Data",
      (Prob_repunit < PofZ_threshold) & (RoSA == "Early") ~ "Mixed Spring",
      (Prob_repunit < PofZ_threshold) & (RoSA == "Late") ~ "Mixed Fall / Late Fall",
      Prob_repunit < PofZ_threshold ~ "Mixed",
      (RoSA == "Early") & (repunit %in% c("fall", "latefall")) ~ "spring",
      (RoSA == "Late") & (repunit == "spring") ~ "Fall",
      (fraction_missing < gsi_missing_threshold) ~ repunit,
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
      (RoSA == "Early") & (repunit %in% c("fall", "latefall")) & (Prob_repunit >= PofZ_threshold) & ((n_miss_loci / (n_miss_loci + n_non_miss_loci)) < gsi_missing_threshold) ~ "Feather River Spring",
      (RoSA == "Late") & (repunit == "spring") & (Prob_repunit > PofZ_threshold) ~ NA_character_,
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
    trib_PofZ = if_else(!(final_call %in% c("Missing Data", "Fall / Late Fall", "Spring / Winter")), round(trib_PofZ, digits = 3), NA_real_),
    final_call = str_to_title(final_call)
  )

mix_results_wide_w_extras <- mix_results_wide %>%
  left_join(species_results, by = "SampleID") %>%
  left_join(heterozygosity_results, by = "SampleID") %>%
  left_join(LFAR_results, by = "SampleID") %>%
  mutate(
    final_call = case_when(
      Species == "non-Chinook" ~ "non-Chinook",
      (LFAR_markers_present == FALSE) & (final_call %in% c("Fall", "Latefall")) ~ "Fall / Late Fall", # If LFAR markers are not present and final call is Fall or Latefall, change final call to Fall / Late Fall
      TRUE ~ final_call
    )
  )
write_tsv(
  mix_results_wide_w_extras,
  file = stringr::str_c(project_name, "_summary.tsv")
)
