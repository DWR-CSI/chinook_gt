#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
samplesheet_file <- args[1]
vcf_file <- args[2]
project_name <- args[3]
num_cpus <- as.numeric(args[4])


# Setup ---------------
library(microhaplot)

haplo.read.tbl <- prepHaplotFiles(
    run.label = project_name,
    sam.path = '.',
    label.path = samplesheet_file,
    vcf.path = vcf_file,
    app.path = '.',
    n.jobs = num_cpus
)