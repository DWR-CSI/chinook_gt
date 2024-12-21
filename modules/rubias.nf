process RUN_RUBIAS {
    tag "Running RUBIAS"
    container "docker.io/bnguyen29/r-rubias:1.0.4"
    label 'process_high'
    publishDir "${params.outdir}/${params.project}/rubias", mode: 'copy'

    input:
    path ots28_report
    path unknowns_numgeno
    path baseline
    val panel
    path unknowns_geno

    output:
    path "*_matchy_pairs.tsv", emit: mpairs
    path "*_full_mix_posteriors.tsv", emit: full_mix_posteriors
    path "*_summary.tsv", emit: summary

    script:
    """
    rubias.R $unknowns_numgeno $baseline ${params.project} ${params.rubias_show_missing_data} ${ots28_report} ${panel} $unknowns_geno
    """
}