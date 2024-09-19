process GEN_MHP_SAMPLE_SHEET {
    tag "Generate MHP samplesheet"
    label 'process_xsmall'
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
    path locindex

    output:
    path "*_genotypes.txt", emit: geno
    path "*_numgenotypes.txt", emit: numgeno
    path "*_missingdata.txt", emit: missing
    path "*_locus_indices.csv", emit: new_index
    path "*_numgenotypes_OTS28.txt", emit: numgeno_OTS28

    script:
    """
    haps2geno.R ${hap_file} ${params.project} ${locindex}
    """
}

process CHECK_FILE_UPDATE {
    tag "Comparing loci indices"
    label 'process_xsmall'
    container 'ubuntu:latest'
    input:
    path new_file
    path old_file

    output:
    stdout

    script:
    """
    if [ -f "${old_file}" ]; then
        if cmp -s "${new_file}" "${old_file}"; then
            echo "Locus index files are identical. No update needed for ${old_file}"
        else
            echo "Input and output locus index files are different. ${old_file} should be updated with content from ${new_file}"
        fi
    else
        echo "${old_file} does not exist. It should be created with content from ${new_file}"
    fi
    """
}