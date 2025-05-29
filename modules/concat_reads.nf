process CONCAT_READS {
    tag "Concatenating reads for ${sample_id}"
    container 'docker.io/nfcore/base:2.1'
    label 'process_small'
    
    publishDir "${params.outdir}/${params.project}/concatenated", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    tuple val(sample_id), path("${sample_id}_all_reads.fastq.gz"), emit: concatenated

    script:
    """
    cat ${read_files} > ${sample_id}_all_reads.fastq.gz
    """
}
