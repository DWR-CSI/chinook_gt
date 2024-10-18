#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
samplesheet_file <- args[1]
vcf_file <- args[2]
project_name <- args[3]
num_cpus <- as.numeric(args[4])
reference_name <- args[5]


# Setup ---------------
library(microhaplot)

haplo.read.tbl <- prepHaplotFiles(
    run.label = paste0(project_name, "_", reference_name),
    sam.path = ".",
    label.path = samplesheet_file,
    vcf.path = vcf_file,
    app.path = ".",
    n.jobs = num_cpus
)
