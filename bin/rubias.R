#!/usr/bin/env Rscript
# Setup --------------
library(tidyverse)
library(data.table)
library(rubias)

# Get command-line arguments passed by Nextflow
args <- commandArgs(trailingOnly = TRUE)
unks <- read_csv(args[1]) %>%
    mutate_if(is.double, as.integer) %>%
    rename_at(vars(ends_with(".1")), ~str_remove(., "\\.1$")) %>%
    rename_at(vars(ends_with(".2")), ~str_replace(., "\\.2$", ".1")) %>%
    mutate(sample_type = "mixture", repunit = NA) %>%
    rename(collection = group, indiv = indiv.ID) %>%
    rename_all(~gsub("-", ".", .)) # dashes not accepted in column names
ref_baseline <- read_csv(args[2]) %>%
  mutate_if(is.double, as.integer) %>%
  rename_all(~gsub("-", ".", .))
project_name <- args[3]
reporting_groups <- as.integer(args[4])
if (reporting_groups != 2 && reporting_groups != 4) {
    stop("reporting_groups must be 2 or 4")
}

show_missing_data <- as.logical(args[5])



# Combine unknowns and reference baseline ----------------
unk_match <- unks %>% 
    select(any_of(names(ref_baseline)))

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
  filter(collection != "Coho")



# Matching --------------------
matchy_pairs <- close_matching_samples(D = chinook_all, 
                                       gen_start_col = 5, 
                                       min_frac_non_miss = 0.85, 
                                       min_frac_matching = 0.94)

matchy_pairs_out <- matchy_pairs %>% 
  arrange(desc(num_non_miss), desc(num_match))

write_tsv(matchy_pairs_out, file = stringr::str_c(project_name, "_matchy_pairs.tsv"))

# Estimate mixtures --------------------

mix_est <- infer_mixture(reference = ref_match,
                         mixture = unk_match,
                         gen_start_col = 5)


### now for exporting the results!
### for the "wide" DFs, you'll need to adjust the number of columns
### to match the number of repunits!!


mix_results <- mix_est$indiv_posteriors %>%
  group_by(indiv, mixture_collection) %>%
  filter(PofZ==max(PofZ)) %>%
  arrange(indiv, collection, repunit, PofZ)

mix_est_export <- mix_results[,1:9]

write.table(
    mix_est_export, 
    file = stringr::str_c(project_name, "_mix_estimates.tsv"),
    sep = "\t", 
    row.names = FALSE
    )


### export all the scaled likelihoods for each repunit

mix_results_long <- mix_est$indiv_posteriors[,c(2:3,5)]

mix_results_wide <- mix_results_long %>%
  spread(repunit, PofZ)


# set the number of cols in mix_est_wide to equal the number of repunits in that DF!
# 4 reporting groups for full panel, 2 reporting groups for reduced panel

if (reporting_groups == 4) {
    mix_report <- merge(
    mix_est_export[,c(1:3,5,9)], 
    mix_results_wide[,c(1:5)], 
    by="indiv") 
    write.table(
        mix_report, 
        file = stringr::str_c(project_name, "_transition_4repunits_mixreport.csv"),
        sep = ",", 
        row.names = FALSE
        )
} else {
    mix_report <- merge(mix_est_export[,c(1:3,5,9)], mix_results_wide[,c(1:3)], by="indiv") #use me for 2 reporting groups
    write.table(
        mix_report, 
        file = stringr::str_c(project_name, "_reduced_panel_2repunits_mixreport.csv"),
        sep = ",", 
        row.names = FALSE
        )
}
