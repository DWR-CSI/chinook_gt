process STRUCTURE_ROSA_REPORT {
    tag "Synthesizing run IDs"
    container 'docker.io/bnguyen29/r-rubias:1.0.4'
    label 'process_high'
    publishDir "${params.outdir}/${params.project}/report", mode: 'copy'

    input:
    path structure_output
    val project
    path structure_input

    output:
    path "*_ots28_report.tsv", emit: ots28_report

    script:
    """
    parse_structure_RoSA.R ${structure_output} ${params.project} ${structure_input} ${params.ots28_missing_threshold}
    """
}