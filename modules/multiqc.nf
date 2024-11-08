process MULTIQC {
    tag "MULTIQC on ${params.project}"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/multiqc:1.23--pyhdfd78af_0':
    'quay.io/biocontainers/multiqc:1.23--pyhdfd78af_0' }"

    publishDir "${params.outdir}/${params.project}/multiqc", mode: 'copy'
    
    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc .
    """
}
