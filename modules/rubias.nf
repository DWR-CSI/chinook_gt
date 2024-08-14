process RUN_RUBIAS {
    tag "Running RUBIAS"
    container "docker.io/bnguyen29/r-rubias:1.0.4"
    label 'process_medium'
    publishDir "${params.outdir}/${params.project}/rubias", mode: 'copy'

    input:
    path unknowns
    path baseline

    output:
    path "*_matchy_pairs.tsv", emit: mpairs
    path "*_mix_estimates.csv", emit: mix_estimates

    script:
    """
    rubias.R $unknowns $baseline ${params.project} ${params.reporting_groups}
    """
}