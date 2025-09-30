#!/usr/bin/env Rscript
library(tidyverse)
library(microhaplotopia)
args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]

# Filtering thresholds
min_hap_depth <- as.numeric(args[2])
min_depth <- as.numeric(args[3])
allele_balance_param <- as.numeric(args[4])

# Custom functions
find_missing_samples2 <- function (raw_data, filtered_data) 
{
  missing_samples <- raw_data %>% filter(!indiv.ID %in% filtered_data$indiv.ID) %>% 
    distinct(indiv.ID, .keep_all = TRUE) %>% dplyr::select(indiv.ID, group)
  missing_samples
}

# 1. Find the correct .rds file
dirFiles <- list.files()
rds.file <- grep(".rds", dirFiles)
pos.rds.file <- grep("_posinfo.rds", dirFiles)
annotate.rds.file <- grep("_annotate.rds", dirFiles)
rds.file <- setdiff(rds.file, c(pos.rds.file, annotate.rds.file))
input_file <- dirFiles[rds.file[1]]  # Select the first matching file

# 2. Load and preprocess the data
mhp_RDS_file <- readRDS(input_file)
hap <- mhp_RDS_file %>%
  rename(
    indiv.ID = id
    
  ) %>% 
  dplyr::filter(haplo != "haplo") %>%
  mutate(
    indiv.ID = str_remove(indiv.ID,"_*$") # Remove any trailing underscores from individual IDs
  )

unfiltered_output_filename <- paste0(project_name, "_observed_unfiltered_haplotype.csv")
write_csv(hap, unfiltered_output_filename) # Don't need to drop first column downstream if using this

hap_fil <- filter_raw_microhap_data( # Filter based on depth and allele balance
  hap,
  haplotype_depth = min_hap_depth,
  total_depth = min_depth, 
  allele_balance = allele_balance_param
)

nxa <- find_NXAlleles(hap_fil) # Find and drop N/X alleles
write_csv(nxa, file = paste0(project_name, "_nxa.csv"))

nxaCount <- nxa %>% group_by(locus) %>% summarise(count = n()) %>% arrange(desc(count)) #modified from Anthony's code to count the number of loci impacted by extra alleles per individual
write_csv(nxaCount, file = paste0(project_name, "_nxa_count.csv"))

hap_fil_nxa <- nxa %>% 
  select(group, indiv.ID, locus) %>% 
  distinct() %>% 
  anti_join(hap_fil, .)

# Find and drop extra aleles
xtralleles <- find_contaminated_samples(hap_fil_nxa)
write_csv(xtralleles, file = paste0(project_name, "_xtraalleles.csv"))

# Extra diagnostic files
indiv <- xtralleles %>% group_by(indiv.ID) %>% summarise(count = n_distinct(locus)) #modified from Anthony's code to count the number of loci impacted by extra alleles per individual
write_csv(indiv, file = paste0(project_name, "_xtraalleles_individuals.csv"))

locus <- xtralleles %>% group_by(locus) %>% summarise(count = n_distinct(indiv.ID)) %>% arrange(desc(count)) #modified from Anthony's code to count the number of individuals per locus that were impacted by extra alleles
write_csv(locus, file = paste0(project_name, "_xtraalleles_locus.csv"))

hap_fil1 <- hap_fil_nxa %>% 
  anti_join(xtralleles)

loc_depth <- summarize_data(
  datafile = hap_fil1,
  group_var = "locus") %>% 
  arrange(., n_samples)
write_csv(loc_depth, file = paste0(project_name, "_locus_depth.csv"))

ind_depth <- summarize_data(
  datafile = hap_fil1,
  group_var = "indiv.ID") %>% 
  arrange(., mean_depth)
write_csv(ind_depth, file = paste0(project_name, "_individual_depth.csv"))

# Drop duplicated samples
duplicate_samples <- find_duplicates(hap_fil1)
if (!is.null(duplicate_samples)){
  hap_fil1 <- resolve_duplicate_samples(hap_fil1, resolve = "drop")
}

# Add second allele for homozygous haplotypes when missing
hap_final <- add_hom_second_allele(hap_fil1)

# Format for RUBIAS/CKMR
haps_2col <- mhap_transform(
  long_genos = hap_final,
  program = "rubias"
)

# append "_1" and "_2" to the locus names
suffs <- c(1,2)
locs <- colnames(haps_2col)[-1]
addnums <- as_tibble(cbind(locname = locs, suffix = suffs)) %>% 
  mutate(twocol = paste(locs, suffs, sep = ".")) %>% 
  pull(twocol)

haps_2col_final <- haps_2col
names(haps_2col_final) <- c("indiv.ID", addnums)

# add back in missing individuals
missing_samples <- find_missing_samples2(hap, hap_fil1)
haps_2col_final <- haps_2col_final %>% add_row(indiv.ID=missing_samples$indiv.ID) %>%
  mutate(group = "ND") %>%
  relocate(group)

# write final genotype file
write_csv(haps_2col_final, file=paste0(project_name, "_filtered_haplotypes.csv"))
