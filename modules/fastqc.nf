process FASTQC {
    tag "FASTQC on $sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    publishDir "${params.outdir}/${params.project}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -q $reads
    """
}