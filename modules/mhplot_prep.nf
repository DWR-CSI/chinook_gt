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
    echo "$aligned_sam_files" >> ${params.project}_mhp_samplesheet.tsv
    """
}