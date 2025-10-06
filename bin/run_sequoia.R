#!/usr/bin/env Rscript
#' Sequoia Parentage Analysis Script
#'
#' Performs parentage analysis using the Sequoia R package.
#' Converts genotype data to major allele counts and runs parentage assignment.
#'
#' Usage:
#'   Rscript run_sequoia.R <sequoia_mode> <ref_genotypes> <ref_lifehistory> \
#'                        <offspring_genotypes> <offspring_birthyear> \
#'                        <offspring_minBY> <offspring_maxBY> <max_age> \
#'                        <project_name> [missing_threshold]
#'
#' Arguments:
#'   sequoia_mode           Mode for Sequoia (e.g., "par" for parentage)
#'   ref_genotypes          Path to reference parental genotypes file
#'   ref_lifehistory        Path to reference parental life history file
#'   offspring_genotypes    Path to offspring genotypes file
#'   offspring_birthyear    Birth year for offspring (integer or "unknown")
#'   offspring_minBY        Minimum birth year for offspring
#'   offspring_maxBY        Maximum birth year for offspring
#'   max_age                Maximum age of offspring
#'   project_name           Project name for output files
#'   missing_threshold      Threshold for missing data filtering (default: 0.5)
#'
#' Outputs:
#'   <project_name>_allele_dictionary.txt   Allele frequency dictionary
#'   <project_name>_sequoia_output.rds       Full Sequoia output object
#'   <project_name>_parentage_results.txt    Parentage assignment results

# Set up libraries
library("tidyverse")
library("sequoia")

# Helper functions ---------------------

#' Check if file exists and is readable
#' @param file_path Path to the file
#' @param description Description of the file for error messages
check_file_exists <- function(file_path, description) {
    if (!file.exists(file_path)) {
        stop("File not found: ", file_path, " (", description, ")")
    }
    if (!file.access(file_path, mode = 4) == 0) {
        stop("File not readable: ", file_path, " (", description, ")")
    }
}

#' Safely convert string to integer with error handling
#' @param value String value to convert
#' @param param_name Parameter name for error messages
#' @param allow_na Whether NA values are allowed
safe_as_integer <- function(value, param_name, allow_na = TRUE) {
    if (is.na(value) || value == "unknown" || value == "") {
        if (allow_na) {
            return(NA_integer_)
        } else {
            stop("Invalid ", param_name, ": cannot be missing or 'unknown'")
        }
    }

    result <- suppressWarnings(as.integer(value))
    if (is.na(result)) {
        stop("Invalid ", param_name, ": '", value, "' is not a valid integer")
    }
    return(result)
}

#' Safely convert string to numeric with error handling
#' @param value String value to convert
#' @param param_name Parameter name for error messages
#' @param default_value Default value if conversion fails
safe_as_numeric <- function(value, param_name, default_value = NULL) {
    result <- suppressWarnings(as.numeric(value))
    if (is.na(result)) {
        if (!is.null(default_value)) {
            return(default_value)
        }
        stop("Invalid ", param_name, ": '", value, "' is not a valid number")
    }
    return(result)
}

# Functions ---------------------
create_allele_dictionary <- function(combined_data) {
    cat("Creating allele dictionary from combined data...\n")

    # Get locus names (remove .1/.2 suffix)
    all_cols <- colnames(combined_data)[colnames(combined_data) != "SAMPLE_ID"]
    locus_names <- unique(gsub("\\.[12]$", "", all_cols))

    cat("Found", length(locus_names), "loci\n")

    # Initialize dictionary using tibble
    dictionary <- tibble(
        locus_name = character(),
        major_allele = character(),
        major_frequency = numeric(),
        total_alleles = integer(),
        missing_rate = numeric()
    )

    # Process each locus
    for (locus in locus_names) {
        col1_name <- paste0(locus, ".1")
        col2_name <- paste0(locus, ".2")

        # Skip if columns don't exist
        if (!col1_name %in% all_cols || !col2_name %in% all_cols) {
            next
        }

        # Get allele data
        col1_data <- as.character(combined_data[[col1_name]])
        col2_data <- as.character(combined_data[[col2_name]])

        # Combine all alleles
        all_alleles <- c(col1_data, col2_data)
        valid_alleles <- all_alleles[!is.na(all_alleles) &
            all_alleles != "*" &
            all_alleles != "ND" &
            all_alleles != ""]

        missing_count <- length(all_alleles) - length(valid_alleles)
        missing_rate <- missing_count / length(all_alleles)

        # Find major allele
        if (length(valid_alleles) > 0) {
            allele_counts <- table(valid_alleles)
            major_allele <- names(which.max(allele_counts))
            major_frequency <- max(allele_counts) / sum(allele_counts)
            total_alleles <- length(unique(valid_alleles))
        } else {
            major_allele <- NA
            major_frequency <- 0
            total_alleles <- 0
        }

        # Add to dictionary
        dictionary <- dictionary %>%
            add_row(
                locus_name = locus,
                major_allele = major_allele,
                major_frequency = major_frequency,
                total_alleles = total_alleles,
                missing_rate = missing_rate
            )
    }

    cat("Dictionary created for", nrow(dictionary), "loci\n")
    cat("Mean major allele frequency:", round(mean(dictionary$major_frequency, na.rm = TRUE), 3), "\n")
    cat("Loci with >20% missing data:", sum(dictionary$missing_rate > 0.2, na.rm = TRUE), "\n")

    return(dictionary)
}

convert_to_major_allele_counts <- function(genotype_data, allele_dictionary) {
    cat("Converting genotype data to major allele counts...\n")

    # Start with SAMPLE_ID
    result <- genotype_data %>% select(SAMPLE_ID)

    # Process each locus in the dictionary
    for (i in 1:nrow(allele_dictionary)) {
        locus_name <- allele_dictionary$locus_name[i]
        major_allele <- allele_dictionary$major_allele[i]

        col1_name <- paste0(locus_name, ".1")
        col2_name <- paste0(locus_name, ".2")

        # Skip if columns don't exist in data
        if (!col1_name %in% colnames(genotype_data) ||
            !col2_name %in% colnames(genotype_data)) {
            next
        }

        # Get allele data
        col1_data <- as.character(genotype_data[[col1_name]])
        col2_data <- as.character(genotype_data[[col2_name]])

        # Count major alleles using vectorized approach
        major_counts <- map2_dbl(col1_data, col2_data, ~ {
            # Check for missing data
            if (is.na(.x) || is.na(.y) ||
                .x %in% c("*", "ND", "") ||
                .y %in% c("*", "ND", "")) {
                return(-9)
            }

            # Check if major allele is available
            if (is.na(major_allele)) {
                return(-9)
            }

            # Count major alleles
            sum(c(.x, .y) == major_allele)
        })

        # Add to result
        result[[locus_name]] <- major_counts
    }

    cat("Converted", nrow(result), "individuals across", ncol(result) - 1, "loci\n")
    return(result)
}


# Parse and setup command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 12) {
    stop("Usage: run_sequoia.R <sequoia_mode> <ref_genotypes> <ref_lifehistory> <offspring_genotypes> <offspring_birthyear> <offspring_minBY> <offspring_maxBY> <max_age> <project_name> [missing_threshold] <species_min_repro_age> <species_max_repro_age>")
}

# Validate and parse arguments
sequoia_mode <- args[1]
ref_genotypes_path <- args[2]
ref_lifehistory_path <- args[3]
offspring_genotypes_path <- args[4]
offspring_birthyear <- safe_as_integer(args[5], "offspring_birthyear")
offspring_minBY <- safe_as_integer(args[6], "offspring_minBY")
offspring_maxBY <- safe_as_integer(args[7], "offspring_maxBY")
max_age <- safe_as_integer(args[8], "max_age", allow_na = FALSE) # max age of offspring
project_name <- args[9]
missing_thresholds <- if (length(args) >= 10) {
    safe_as_numeric(args[10], "missing_threshold", 0.5)
} else {
    0.5
}
species_min_repro_age <- args[11] %>%
    as.integer()
species_max_repro_age <- args[12] %>%
    as.integer()

# Validate reproductive age parameters
if (is.na(species_min_repro_age) || is.na(species_max_repro_age)) {
    stop("species_min_repro_age and species_max_repro_age must be valid integers")
}
if (species_min_repro_age >= species_max_repro_age) {
    stop("species_min_repro_age must be less than species_max_repro_age")
}
if (species_min_repro_age < 0 || species_max_repro_age < 0) {
    stop("Reproductive ages must be non-negative")
}

# Validate file existence
check_file_exists(ref_genotypes_path, "reference parental genotypes")
check_file_exists(ref_lifehistory_path, "reference parental life history")
check_file_exists(offspring_genotypes_path, "offspring genotypes")

# Read data files
reference_parental_genotypes <- ref_genotypes_path %>%
    read_tsv(col_types = cols(
        .default = col_character(),
        SAMPLE_ID = col_character()
    ))

reference_parental_lifehistory <- ref_lifehistory_path %>%
    read_tsv(col_types = cols(
        .default = col_character(),
        ID = col_character(),
        Sex = col_integer(),
        BirthYear = col_integer(),
        BY.min = col_integer(),
        BY.max = col_integer(),
        Year.last = col_integer()
    ))

offspring_genotypes <- offspring_genotypes_path %>%
    read_csv(col_types = cols(
        .default = col_character(),
        indiv.ID = col_character()
    )) %>%
    {
        if ("group" %in% colnames(.)) select(., -group) else .
    } %>%
    rename(SAMPLE_ID = indiv.ID)

# Main execution
offspring_lh <- tibble(
    ID = offspring_genotypes$SAMPLE_ID,
    Sex = as.integer(3), # Unknown
    BirthYear = offspring_birthyear,
    BY.min = offspring_minBY,
    BY.max = offspring_maxBY,
    Year.last = NA_integer_
) %>%
    mutate(
        BY.min = case_when(
            is.na(BirthYear) & is.na(BY.min) & !is.na(BY.max) ~ BY.max - max_age, # Default to 8 years before maxBY if no birthyear or minBY is provided, Chinook salmon typical max lifespan.
            TRUE ~ BY.min
        )
    )

combined_genotypes <- bind_rows(
    reference_parental_genotypes,
    offspring_genotypes
)

# Get loci removal regex from environment variable
loci_removal_regex <- Sys.getenv("LOCI_REMOVAL_REGEX")

# Handle empty regex - if empty, set to pattern that matches nothing
if (loci_removal_regex == "") {
    loci_removal_regex <- "^$"
    cat("No loci removal regex specified - no loci will be removed\n")
}

# Print matching columns to be removed
cols_to_remove <- tryCatch(
    names(combined_genotypes)[grepl(loci_removal_regex,
                                     names(combined_genotypes),
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

# Remove columns using pre-calculated list
if (length(cols_to_remove) > 0) {
    combined_genotypes <- combined_genotypes %>%
        select(-all_of(cols_to_remove))
}

allele_dict <- create_allele_dictionary(combined_genotypes)
write_tsv(allele_dict, paste0(project_name, "_allele_dictionary.txt"))
major_allele_counts <- convert_to_major_allele_counts(combined_genotypes, allele_dict)

geno_matrix <- major_allele_counts %>%
    select(-SAMPLE_ID) %>%
    as.matrix()
rownames(geno_matrix) <- major_allele_counts$SAMPLE_ID
write_tsv(
    as.data.frame(major_allele_counts),
    paste0(project_name, "_sequoia_genotype_matrix.txt"),
    col_names = TRUE,
)



combined_lh <- bind_rows(
    reference_parental_lifehistory,
    offspring_lh
)
write_tsv(
    combined_lh,
    paste0(project_name, "_sequoia_lifehistory_data.txt"),
    col_names = TRUE
)

# Remove individuals with missing data above threshold
cat("Filtering individuals with >", missing_thresholds * 100, "% missing data...\n")
missing_rates <- rowMeans(is.na(geno_matrix))
filtered_geno_matrix <- geno_matrix[missing_rates < missing_thresholds, ]
cat("Retained", nrow(filtered_geno_matrix), "of", nrow(geno_matrix), "individuals after filtering\n")

# Individuals with both life history and genotype data
retained_ids <- intersect(
    rownames(filtered_geno_matrix),
    combined_lh$ID
)

final_geno_matrix <- filtered_geno_matrix %>%
    .[rownames(.) %in% retained_ids, ]
final_lh <- combined_lh %>%
    filter(ID %in% retained_ids) %>%
    as.data.frame() # Cannot be a tibble or else you will get an error from Sequoia
if (nrow(final_geno_matrix) == 0 || nrow(final_lh) == 0) {
    stop("No individuals with both genotype and life history data after filtering.")
}
if (nrow(final_geno_matrix) != nrow(final_lh)) {
    stop("Mismatch between genotype and life history data after filtering. Check data.")
}
# Run Sequoia
SeqOUT <- sequoia(
    GenoM = final_geno_matrix,
    LifeHistData = final_lh,
    Err = 0.005,
    Module = sequoia_mode,
    args.AP = list(
        Discrete = TRUE,
        MinAgeParent = species_min_repro_age, # Minimum age of parents
        MaxAgeParent = species_max_repro_age # Maximum age of parents
    ),
    quiet = "verbose",
    Plot = FALSE
)

saveRDS(SeqOUT, paste0(project_name, "_sequoia_output.rds"))

parentage_results <- SeqOUT$Pedigree

if (!is.null(parentage_results)) {
    write_tsv(parentage_results, paste0(project_name, "_parentage_results.txt"))
} else {
    cat("No pedigree results - likely running in preprocessing mode only\n")
    # Create empty parentage results file to satisfy output requirements
    empty_parentage <- data.frame(
        ID = character(0),
        Dam = character(0),
        Sire = character(0),
        stringsAsFactors = FALSE
    )
    write_tsv(empty_parentage, paste0(project_name, "_parentage_results.txt"))
}
