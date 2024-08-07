#!/usr/bin/env Rscript
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]

# 1. Find the correct .rds file
dirFiles <- list.files()
rds.file <- grep(".rds", dirFiles)
pos.rds.file <- grep("_posinfo.rds", dirFiles)
annotate.rds.file <- grep("_annotate.rds", dirFiles)
rds.file <- setdiff(rds.file, c(pos.rds.file, annotate.rds.file))
input_file <- dirFiles[rds.file[1]]  # Select the first matching file

# 2. Load and preprocess the data
table.out <- readRDS(input_file) %>%
ungroup() %>%
mutate(id = as.character(id),
        locus = as.character(locus),
        group = as.character(group))

if ("sum.Phred.C" %in% colnames(table.out)) {
table.out <- table.out %>% select(-sum.Phred.C, -max.Phred.C)
}

# 3. Apply minimal filtering (in this case, we're not filtering at all to get all unfiltered data)
haplo.all <- table.out

# 4. Rename the id column and write to CSV
haplo.all <- haplo.all %>% rename("indiv.ID" = id)
output_filename <- paste0(project_name, "_observed_unfiltered_haplotype.csv")
write.csv(haplo.all, output_filename) # Don't need to drop first column downstream if using this