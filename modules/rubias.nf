process RUN_RUBIAS {
    tag "Genetic stock assignment: RUBIAS"
    container "docker.io/bnguyen29/r-rubias:1.0.4"
    label 'process_high'
    publishDir "${params.outdir}/${params.project}/rubias", mode: 'copy'

    input:
    path ots28_report
    path baseline
    val panel
    path unknowns_geno

    output:
    path "*_matchy_pairs.tsv", emit: mpairs
    path "*_full_mix_posteriors.tsv", emit: full_mix_posteriors
    path "*_summary.tsv", emit: summary

    script:
    """
    rubias.R $baseline ${params.project} ${params.rubias_show_missing_data} ${ots28_report} ${panel} $unknowns_geno ${params.ots28_missing_threshold} ${params.gsi_missing_threshold} ${params.pofz_threshold}
    """
}