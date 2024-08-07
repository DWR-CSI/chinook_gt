process GEN_MHP_SAMPLE_SHEET {
    tag "Generate MHP samplesheet"
    label 'process_small'
    container 'ubuntu:latest'

    publishDir "${params.outdir}/${params.project}", mode: 'copy', pattern: '*_mhp_samplesheet.tsv'

    input:
    path aligned_sam_files

    output:
    path "${params.project}_mhp_samplesheet.tsv", emit: mhp_samplesheet

    script:
    """
    #!/bin/bash
    for file in \$(ls \$aligned_sam_files); do
        echo -e "\$file\t\${file::-8}\tNA" >> "${params.project}_mhp_samplesheet.tsv"
    done
    sort -u -o "${params.project}_mhp_samplesheet.tsv" "${params.project}_mhp_samplesheet.tsv"
    """
}

process PREP_MHP_RDS {
    tag "Prepare Microhaplotype RDS files"
    label 'process_high'
    container 'docker.io/bnguyen29/r-rubias:1.0.4'
    publishDir "${params.outdir}/${params.project}/mhp_rds", mode: 'copy'

    input:
    path samplesheet
    path vcf_file
    path sam_files

    output:
    path "*.rds", emit: rds

    script:
    """
    prep_mhp_rds.R ${samplesheet} ${vcf_file} ${params.project} ${task.cpus}
    """
}

process GEN_HAPS {
    tag "Generate haplotypes"
    label 'process_small'
    container 'docker.io/rocker/tidyverse:4.4.1'
    publishDir "${params.outdir}/${params.project}/haplotypes", mode: 'copy'

    input:
    path rds_files

    output:
    path "*_observed_unfiltered_haplotype.csv", emit: haps

    script:
    """
    gen_haps.R ${params.project}
    """
}

process HAP2GENO {
    tag "Convert haplotypes to genotypes"
    label 'process_small'
    container 'docker.io/rocker/tidyverse:4.4.1'
    publishDir "${params.outdir}/${params.project}/genotypes", mode: 'copy'

    input:
    path hap_file

    output:
    path "*_genotypes.txt", emit: geno
    path "*_numgenotypes.txt", emit: numgeno
    path "*_missingdata.txt", emit: missing

    script:
    """
    haps2geno.R ${hap_file} ${params.project}
    """
}