#!/usr/bin/env Rscript
library(tidyverse)
library(microhaplotopia)
args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]

# Filtering thresholds
min_hap_depth <- as.numeric(args[2])
min_depth <- as.numeric(args[3])
allele_balance_param <- as.numeric(args[4])

# 1. Find the correct .rds file
dirFiles <- list.files()
rds.file <- grep(".rds", dirFiles)
pos.rds.file <- grep("_posinfo.rds", dirFiles)
annotate.rds.file <- grep("_annotate.rds", dirFiles)
rds.file <- setdiff(rds.file, c(pos.rds.file, annotate.rds.file))
input_file <- dirFiles[rds.file[1]]  # Select the first matching file

# 2. Load and preprocess the data
mhp_RDS_file <- readRDS(input_file) %>%
hap <- mhp_RDS_file %>%
  rename(
    indiv.ID = id
    
  ) %>% 
  dplyr::filter(haplo != "haplo")

hap_fil <- filter_raw_microhap_data( # Filter based on depth and allele balance
  hap,
  haplotype_depth = hapDepth,
  total_depth = totDepth, 
  allele_balance = alleleBalance
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
write_csv(xtralleles, file = paste0(project_name, "_xtralleles.csv"))

# Extra diagnostic files
indiv <- xtralleles %>% group_by(indiv.ID) %>% summarise(count = n_distinct(locus)) #modified from Anthony's code to count the number of loci impacted by extra alleles per individual
write_csv(indiv, file = paste0(project_name, "_xtralleles_individuals.csv"))

locus <- xtralleles %>% group_by(locus) %>% summarise(count = n_distinct(indiv.ID)) %>% arrange(desc(count)) #modified from Anthony's code to count the number of individuals per locus that were impacted by extra alleles
write_csv(locus, file = paste0(project_name, "_xtralleles_locus.csv"))

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

haplo.all <- haplo.all %>% rename("indiv.ID" = id)
output_filename <- paste0(project_name, "_observed_unfiltered_haplotype.csv")
write.csv(haplo.all, output_filename) # Don't need to drop first column downstream if using this

# 5. Filtered haplotype production
