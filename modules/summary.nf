process RUBIAS_SUMMARY {
    tag "Generating summary statistics for RUBIAS: $reference"
    label 'process_small'
    container 'docker.io/rocker/tidyverse:4.5.2'
    publishDir "${params.outdir}/${params.project}/summary", mode: 'copy'

    input:
    path rubias_results

    output:
    path "*_rubias_summary.tsv", emit: rubias_summary

    script:
    """
    summarize_results.R rubias ${params.project} ${rubias_results}
    """
}

process CKMRSIM_RUBIAS_SUMMARY {
    tag "Generating summary statistics for CKMRsim RUBIAS: $reference"
    label 'process_small'
    container 'docker.io/rocker/tidyverse:4.5.2'
    publishDir "${params.outdir}/${params.project}/summary", mode: 'copy'

    input:
    path rubias_results
    path ckmrsim_results

    output:
    path "*_ckmrsim_rubias_summary.tsv", emit: ckmrsim_rubias_summary

    script:
    """
    summarize_results.R ckmrsim_rubias ${params.project} ${rubias_results} ${ckmrsim_results}
    """
}