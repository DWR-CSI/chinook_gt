process RUN_ID {
    tag "Synthesizing run IDs"
    container 'docker.io/bnguyen29/r-rubias:1.0.4'
    label 'process_high'
    publishDir "${params.outdir}/${params.project}/report", mode: 'copy'

    input:
    path structure_output
    path gsi_output
    path structure_input

    output:
    path "*_run_id.txt", emit: run_id
    path "*_ots28_report.tsv", emit: ots28_report

    script:
    """
    run_id.R ${structure_output} ${gsi_output} ${params.project} ${structure_input}
    """
}