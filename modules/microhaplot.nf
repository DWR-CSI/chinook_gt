process GEN_MHP_SAMPLE_SHEET {
    tag "Generate MHP samplesheet: $reference"
    label 'process_xsmall'
    container 'docker.io/nfcore/base:2.1'

    publishDir "${params.outdir}/${params.project}/mhp/${reference}", mode: 'copy', pattern: '*_mhp_samplesheet.tsv'

    input:
    tuple val(sample_id), val(reference), path(vcf), path(aligned_sam_files)

    output:
    tuple val(reference), path("${params.project}_${reference}_mhp_samplesheet.tsv"), path(vcf), path(aligned_sam_files), emit: mhp_samplesheet

    script:
    """
    #!/bin/bash
    for sam_file in ${aligned_sam_files}; do
        if [[ "\$sam_file" == *.sam ]]; then
            base_name=\$(basename "\$sam_file" ${reference}_aln.sam)
            echo -e "\$sam_file\t\$base_name\tNA" >> "${params.project}_${reference}_mhp_samplesheet.tsv"
        fi
    done
    sort -u -o "${params.project}_${reference}_mhp_samplesheet.tsv" "${params.project}_${reference}_mhp_samplesheet.tsv"
    """
}

process PREP_MHP_RDS {
    tag "Prepare Microhaplotype RDS files: $reference"
    label 'process_high'
    container 'docker.io/bnguyen29/r-rubias:1.0.4'
    publishDir "${params.outdir}/${params.project}/mhp/${reference}", mode: 'copy'

    input:
    tuple val(reference), path(samplesheet), path(vcf_file), path(sam_files)

    output:
    tuple val(reference), path("*.rds"), emit: rds

    script:
    """
    # Needs to take in reference info and output reference info in the RDS file
    prep_mhp_rds.R ${samplesheet} ${vcf_file} ${params.project} ${task.cpus} ${reference}
    """
}

process GEN_HAPS {
    tag "Generate haplotypes: $reference"
    label 'process_small'
    container 'docker.io/rocker/tidyverse:4.4.1'
    publishDir "${params.outdir}/${params.project}/haplotypes", mode: 'copy'

    input:
    tuple val(reference), path(rds_files)

    output:
    tuple val(reference), path("*_observed_unfiltered_haplotype.csv"), emit: haps

    script:
    """
    gen_haps.R ${params.project}_${reference}
    """
}

process HAP2GENO {
    tag "Convert haplotypes to genotypes: $reference"
    label 'process_small'
    container 'docker.io/rocker/tidyverse:4.4.1'
    publishDir "${params.outdir}/${params.project}/genotypes", mode: 'copy'

    input:
    tuple val(reference), path(hap_file)
    path locindex

    output:
    path "*_genotypes.txt", emit: geno
    path "*_numgenotypes.txt", emit: numgeno
    path "*_missingdata.txt", emit: missing
    path "*_locus_indices.csv", emit: new_index
    path "*_numgenotypes_OTS28.txt", emit: numgeno_OTS28

    script:
    """
    haps2geno.R ${hap_file} ${params.project}_${reference} ${locindex}
    """
}

process CHECK_FILE_UPDATE {
    tag "Comparing loci indices"
    label 'process_xsmall'
    container 'docker.io/nfcore/base:2.1'
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