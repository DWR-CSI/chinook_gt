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