#!/usr/bin/env Rscript
library(tidyverse)
library(CKMRsim)

# Load arguments
args <- commandArgs(trailingOnly = TRUE)
unknown_genotypes_raw <- args[1] %>%
  read_csv(col_types = cols(.default = col_character())) %>%
  select(-'group') %>%
  rename(SAMPLE_ID = 'indiv.ID')
parents_genotypes_raw <- args[2] %>%
  read_csv(col_types = cols(.default = col_character()))
logl_threshold <- args[3]  # Check if it is numeric or a character string. If it is 'auto', we will use automatic thresholding.
if (str_to_lower(logl_threshold) != "auto") {
  logl_threshold <- as.numeric(logl_threshold)
} else {
  logl_threshold <- "auto" # not yet supported
  stop("Automatic log-likelihood thresholding not yet supported in this script.")
}
project_name <- args[4]
# allele_freqs input file should be formatted with columns: Chrom, Locus, Pos, Allele, LocIdx, AlleIdx, Freq
min_loci_threshold <- args[5] %>% as.integer() # Need to add default param setting to main.nf, nextflow schema, and example params.
extra_genos_long <- args[6] %>% # extra long-format genotype data for calculating allele frequencies
    readRDS()

# Get loci removal regex from environment variable
loci_removal_regex <- Sys.getenv("LOCI_REMOVAL_REGEX")
# Handle empty regex - if empty, set to pattern that matches nothing
if (loci_removal_regex == "") {
    loci_removal_regex <- "^$"
    cat("No loci removal regex specified - no loci will be removed\n")
}

# Functions
index_markers <- function(M) {
  tmp <- M %>%
    dplyr::distinct(Chrom, Locus) %>%
    dplyr::count(Locus) %>%
    dplyr::filter(n > 1)
  if (nrow(tmp) > 0) {
    dupies <- tmp$Locus
    stop("Locus names must be globally unique.  These are not: ", paste(dupies, collapse = ", "))
  }
  M %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Chrom, Pos, desc(Freq)) %>%
    dplyr::mutate(locidx = as.integer(factor(Locus, levels = unique(Locus)))) %>%
    dplyr::group_by(Chrom, Locus) %>%
    dplyr::mutate(
      alleidx = as.integer(factor(Allele, levels = unique(Allele))),
      newfreq = Freq/sum(Freq)
      ) %>%
    dplyr::select(-Freq) %>%
    rename(Freq = newfreq, AlleIdx = alleidx, LocIdx = locidx) %>%
    dplyr::ungroup()
}

reshape_paired_genotypes <- function(geno_df) {

  # Extract sample IDs
  sample_ids <- geno_df$SAMPLE_ID

  # Get column names (excluding SAMPLE_ID)
  geno_cols <- colnames(geno_df)[-1]

  # Extract locus names (remove .1 or .2 suffix)
  locus_names <- unique(gsub("\\.[12]$", "", geno_cols))

  # Create long format data frame
  long_list <- list()

  for (i in seq_along(sample_ids)) {
    indiv_id <- sample_ids[i]

    for (locus in locus_names) {
      allele1_col <- paste0(locus, ".1")
      allele2_col <- paste0(locus, ".2")

      allele1 <- as.character(geno_df[[allele1_col]][i])
      allele2 <- as.character(geno_df[[allele2_col]][i])

      # Store genotype
      long_list[[length(long_list) + 1]] <- tibble(
        Indiv = indiv_id,
        Locus = locus,
        gene_copy = 1,
        Allele = allele1
      )

      long_list[[length(long_list) + 1]] <- tibble(
        Indiv = indiv_id,
        Locus = locus,
        gene_copy = 2,
        Allele = allele2
      )
    }
  }

  bind_rows(long_list)
}



# Process genotypes
combined_genotypes_raw <- bind_rows(
  unknown_genotypes_raw,
  parents_genotypes_raw
)
# Remove loci based on regex
cols_to_remove <- tryCatch(
    names(combined_genotypes_raw)[grepl(loci_removal_regex,
                                     names(combined_genotypes_raw),
                                     perl = TRUE)],
    error = function(e) {
        cat("Invalid loci removal regex pattern: ", loci_removal_regex, "\n")
        cat("Error: ", e$message, "\n")
        cat("No loci will be removed\n")
        loci_removal_regex <<- "^$"
        character(0)
    }
)

if (length(cols_to_remove) > 0) {
    cat("Columns to be removed from combined_genotypes:\n")
    cat(paste(cols_to_remove, collapse = ", "), "\n")
}
if (length(cols_to_remove) > 0) {
    combined_genotypes_raw <- combined_genotypes_raw %>%
        select(-all_of(cols_to_remove))
}

offspring_ids <- unknown_genotypes_raw$SAMPLE_ID
parent_ids <- parents_genotypes_raw$SAMPLE_ID

combined_genotypes_long <- reshape_paired_genotypes(combined_genotypes_raw)
# Load and calculate allele frequencies
cat("Loading allele frequencies...\n")

total_genos_long <- bind_rows(combined_genotypes_long, extra_genos_long)
allele_freqs <- total_genos_long %>%
  # Remove missing data (if coded as NA, "", or specific missing code)
  filter(!is.na(Allele), Allele != "", Allele != "0") %>%
  # Count alleles by locus
  group_by(Locus, Allele) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Calculate frequencies within each locus
  group_by(Locus) %>%
  mutate(
    total = sum(count),
    Freq = count / total
  ) %>%
  ungroup() %>%
  select(Locus, Allele, Freq)

# Check allele frequencies
allele_freqs %>%
  group_by(Locus) %>%
  summarise(
    n_alleles = n(),
    freq_sum = sum(Freq)
  ) %>%
  print(n = 20)

# Identify and remove monomorphic loci (only 1 allele)
monomorphic_loci <- allele_freqs %>%
  group_by(Locus) %>%
  summarise(n_alleles = n()) %>%
  filter(n_alleles == 1) %>%
  pull(Locus)

cat("Removing", length(monomorphic_loci), "monomorphic loci\n")

# Filter out monomorphic loci
allele_freqs <- allele_freqs %>%
  filter(!Locus %in% monomorphic_loci)

# Create the long_markers format required by CKMRsim
# This format needs: Chrom, Locus, Pos, Allele, Freq

# For microhaplotypes without known physical positions,
# assign arbitrary chromosome and position values
long_markers_data <- allele_freqs %>%
  mutate(
    Chrom = 1,  # Arbitrary chromosome (use 1 for unlinked analysis)
    Pos = as.numeric(factor(Locus)),  # Arbitrary position based on locus order
  ) %>%
  select(Chrom, Locus, Pos, Allele, Freq)

# Reindex markers (required preprocessing step)
long_markers_data <- index_markers(long_markers_data)

# Check the result
head(long_markers_data, 20)

# Summary statistics
cat("Number of loci:", n_distinct(long_markers_data$Locus), "\n")
cat("Number of total alleles:", nrow(long_markers_data), "\n")
cat("Alleles per locus (range):",
    range(table(long_markers_data$Locus)), "\n")

PO_ckmr <- create_ckmr(
  D = long_markers_data,
  ge_mod_assumed = ge_model_microhap1,
  ge_mod_true = ge_model_microhap1,
  ge_mod_assumed_pars_list = list(miscall_rate = 0.005, dropout_rate = 0.005),
  ge_mod_true_pars_list = list(miscall_rate = 0.005, dropout_rate = 0.005)
)

cat("Number of markers in allele frequencies file:", n_distinct(allele_freqs$Locus), "\n")
cat("Number of total alleles:", nrow(allele_freqs), "\n")
cat("Min and max of alleles per locus (range): ", paste(range(table(allele_freqs$Locus)), collapse = " - "), "\n")


# Remove geno loci that are not in allele frequency marker file
geno_loci <- unique(combined_genotypes_long$Locus)
reference_set_loci <- unique(allele_freqs$Locus)
loci_to_remove <- setdiff(geno_loci, reference_set_loci)
if (length(loci_to_remove) > 0) {
  cat("Removing", length(loci_to_remove), "loci from genotypes that are not in allele frequency reference set.\n")
  cat("Loci removed:", paste(loci_to_remove, collapse = ", "), "\n")
  combined_genotypes_long <- combined_genotypes_long %>%
    filter(!Locus %in% loci_to_remove)
}

# Check for duplicates
duplicates <- combined_genotypes_long %>%
  group_by(Indiv, Locus, gene_copy) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

if(nrow(duplicates) > 0) {
  cat("WARNING: Found", nrow(duplicates), "duplicate entries. Displaying first 20:\n")
  print(head(duplicates, 20))
  cat("WARNING: Removing duplicates by keeping first occurrence. This may indicate data quality issues.\n")
  cat("WARNING: Consider investigating the source of duplicate genotype calls.\n")
  combined_genotypes_long <- combined_genotypes_long %>%
    distinct(Indiv, Locus, gene_copy, .keep_all = TRUE)
}

# Split back into parent and offspring genotype dataframes
parent_genos_long <- combined_genotypes_long %>%
  filter(Indiv %in% parent_ids)
offspring_genos_long <- combined_genotypes_long %>%
  filter(Indiv %in% offspring_ids)
write_tsv(parent_genos_long,
          file = paste0(project_name, "_parent_genotypes_long.tsv"))
write_tsv(offspring_genos_long,
          file = paste0(project_name, "_offspring_genotypes_long.tsv"))
cat("Number of offspring genotypes:", length(unique(offspring_genos_long$Indiv)), "\n")
cat("Number of parent genotypes:", length(unique(parent_genos_long$Indiv)), "\n")

po_results <- pairwise_kin_logl_ratios(
    D1 = parent_genos_long,
    D2 = offspring_genos_long,
    CK = PO_ckmr,
    numer = "PO",
    denom = "U"
)

if (logl_threshold != "auto" && !is.na(as.numeric(logl_threshold))) {
  logl_threshold <- as.numeric(logl_threshold)
  cat("Applying log-likelihood ratio threshold of", logl_threshold, "\n")
  po_results_filtered <- po_results %>%
    filter(logl_ratio >= logl_threshold) %>%
    filter(num_loc >= min_loci_threshold) %>%
    arrange(desc(logl_ratio))
} else if (logl_threshold == "auto") {
    stop("Automatic log-likelihood thresholding not yet supported in this script.")
} else {
  po_results_filtered <- po_results %>%
    filter(num_loc >= min_loci_threshold) %>%
    arrange(desc(logl_ratio))
}

write_tsv(po_results_filtered, file = paste0(project_name, "_PO_results.tsv"))
