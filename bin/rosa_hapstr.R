#!/usr/bin/env Rscript
library(tidyverse)
library(vcfR)

# Functions ---------

# Functions ---------

#' Calculate run timing scores and make genotype calls
#' @param ptypes Dataframe containing haplotype strings and counts
#' @param threshold_het Threshold for calling heterozygous (default 0.5)
#' @param min_markers Minimum number of scored markers required (default 3 and 2 for EL and W respectively)
#' @return Dataframe with run timing calls and confidence scores
analyze_run_timing <- function(ptypes, threshold_het = 0.5, min_EL_markers = 3, min_W_markers = 2) {
    ptypes %>%
        mutate(
            # Calculate total markers scored
            total_scored = nchar(hapstr) - Q_count - WQ_count,
            EL_total_scored = E_count + L_count + H_count,
            Wdiag_total_scored = W_count + N_count + M_count,

            # Calculate proportions for each type
            early_prop = E_count / EL_total_scored,
            late_prop = L_count / EL_total_scored,
            het_prop = H_count / EL_total_scored,
            winter_prop = W_count / Wdiag_total_scored,
            ots28_missing = round(Q_count / (EL_total_scored + Q_count) * 100, digits = 0), # Does not include Winter diag loci

            # Calculate run timing scores
            timing_score = (E_count - L_count) / EL_total_scored,

            # Make run timing calls
            timing_call = case_when(
                EL_total_scored < min_EL_markers ~ "Insufficient_Data",
                het_prop >= threshold_het ~ "Intermediate", # Not distinguished from mixed non-hets
                timing_score > 0.5 ~ "Early",
                timing_score < -0.5 ~ "Late",
                TRUE ~ "Intermediate"
            ),

            # Make winter run calls
            winter_call = case_when(
                Wdiag_total_scored < min_W_markers ~ "Insufficient_Data",
                winter_prop >= 0.5 ~ "Winter",
                winter_prop < 0.2 ~ "Non_Winter",
                TRUE ~ "Possibly_Winter"
            )

            # Maybe calculate timing confidence and winter confidence metrics in the future
        )
}

#' Summarize run timing results for a population
#' @param analyzed_data Output from analyze_run_timing()
#' @return List containing summary statistics and plots
summarize_population <- function(analyzed_data) {
    # Calculate population-level statistics
    summary_stats <- analyzed_data %>%
        summarise(
            n_samples = n(),
            n_early = sum(timing_call == "Early"),
            n_late = sum(timing_call == "Late"),
            n_intermed = sum(timing_call == "Intermediate"),
            n_winter = sum(winter_call == "Winter")
        )

    # Create timing score distribution plot
    timing_plot <- ggplot(analyzed_data, aes(x = timing_score)) +
        geom_histogram(binwidth = 0.1) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red") +
        labs(
            title = "Distribution of Run Timing Scores",
            x = "Timing Score (negative = Late, positive = Early)",
            y = "Count"
        )

    # Create winter proportion plot
    winter_plot <- ggplot(analyzed_data, aes(x = winter_prop)) +
        geom_histogram(binwidth = 0.1) +
        geom_vline(xintercept = c(0.2, 0.5), linetype = "dashed", color = "red") +
        labs(
            title = "Distribution of Winter Run Proportions",
            x = "Proportion of Winter Run Markers",
            y = "Count"
        )

    list(
        summary = summary_stats,
        timing_plot = timing_plot,
        winter_plot = winter_plot
    )
}

# Read in command line arguments ---------
args <- commandArgs(trailingOnly = TRUE)
project_name <- as.character(args[1])
vcf_path <- args[2]
allele_key_path <- args[3]

# Main code ---------
locs <- read_tsv(allele_key_path) %>%
    mutate(
        pheno = case_when(
            (str_detect(CHROMPOS, "winter-greb1l-diagnostic") & pheno == "?") ~ "/",
            (str_detect(CHROMPOS, "winter-greb1l-diagnostic") & pheno == "H") ~ "M",
            TRUE ~ pheno
        )
    )
vcf <- read.vcfR(vcf_path)
vcf_tidy <- vcfR2tidy(vcf)

# Create a complete set of expected CHROMPOS values from locs
expected_chrompos <- locs %>%
    filter(gt_GT_alleles != "REF" & gt_GT_alleles != ".") %>%
    distinct(CHROMPOS)

# Process genotypes with complete locus representation
gtypes <- vcf_tidy$gt %>%
    select(ChromKey, POS, Indiv, gt_GT_alleles, gt_DP, gt_AD) %>%
    left_join(
        vcf_tidy$fix %>%
            select(
                ChromKey, CHROM, POS
            )
    ) %>%
    unite("CHROMPOS", c("CHROM", "POS"), sep = "_", remove = T) %>%
    # First ensure we have all combinations of Indiv and CHROMPOS
    full_join(
        crossing(
            CHROMPOS = expected_chrompos$CHROMPOS,
            Indiv = unique(vcf_tidy$gt$Indiv)
        ),
        by = c("CHROMPOS", "Indiv")
    ) %>%
    filter(CHROMPOS %in% locs$CHROMPOS) %>%
    group_by(Indiv) %>%
    arrange(match(CHROMPOS, locs$CHROMPOS), .by_group = TRUE) %>%
    ungroup()

ptypes <- gtypes %>%
    select(CHROMPOS, Indiv, gt_GT_alleles) %>%
    # No longer need complete() since we did full_join above
    mutate(gt_GT_alleles = if_else(is.na(gt_GT_alleles), ".", gt_GT_alleles)) %>%
    left_join(locs) %>%
    group_by(Indiv) %>%
    distinct(Indiv, CHROMPOS, .keep_all = TRUE) %>%
    arrange(haporder, .by_group = TRUE) %>%
    summarise(hapstr = paste(pheno, collapse = "")) %>%
    mutate(
        # Replace any NA values that might have slipped through
        hapstr = str_replace_all(hapstr, "NA", "?"),
        # Count various allele types
        E_count = str_count(hapstr, "E"), # Early
        L_count = str_count(hapstr, "L"), # Late
        Q_count = str_count(hapstr, "\\?"), # Missing, normal RoSA
        H_count = str_count(hapstr, "H"), # Heterozygous
        W_count = str_count(hapstr, "W"), # Winter configuration
        N_count = str_count(hapstr, "N"), # Non-winter configuration
        M_count = str_count(hapstr, "M"), # M = Mixed = Winter diag loci het
        WQ_count = str_count(hapstr, "/") # Winter diag loci, missing
    )

ptypes_analyzed <- analyze_run_timing(ptypes)

write_tsv(ptypes_analyzed, file = str_c(project_name, "_RoSA_hapstrs.tsv"))

ptypes_report <- ptypes_analyzed %>%
    select(indiv = Indiv, RoSA = timing_call, ots28_missing, hapstr) %>%
    mutate(
        RoSA = case_when(
            RoSA == "Insufficient_Data" ~ "Uncertain",
            TRUE ~ RoSA
        )
    )

write_tsv(ptypes_report, file = str_c(project_name, "_hapstr_ots28_report.tsv"))

# summaries <- summarize_population(ptypes_analyzed)
