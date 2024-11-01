process RUN_RUBIAS {
    tag "Running RUBIAS"
    container "docker.io/bnguyen29/r-rubias:1.0.4"
    label 'process_high'
    publishDir "${params.outdir}/${params.project}/rubias", mode: 'copy'

    input:
    path ots28_report
    path unknowns
    path baseline
    val panel

    output:
    path "*_matchy_pairs.tsv", emit: mpairs
    path "*_full_mix_estimates.tsv", emit: full_mix_estimates
    path "*_summary.tsv", emit: summary

    script:
    """
    rubias.R $unknowns $baseline ${params.project} ${params.reporting_groups} ${params.rubias_show_missing_data} ${ots28_report} ${panel}
    """
}