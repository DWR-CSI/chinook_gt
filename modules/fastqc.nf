process FASTQC {
    tag "FASTQC"
    label 'process_medium'
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    publishDir "${params.outdir}/${params.project}/fastqc", mode: 'copy'

    input:
    path(fastq_files)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    def threads = task.cpus > 1 ? task.cpus : 1
    """
    fastqc -q -t ${threads} $fastq_files

    """
}