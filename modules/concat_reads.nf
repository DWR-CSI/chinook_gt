process CONCAT_READS {
    tag "Concatenating reads for ${sample_id}"
    label 'process_small'
    
    publishDir "${params.outdir}/${params.project}/concatenated", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    tuple val(sample_id), path("${sample_id}_all_reads.fastq.gz"), emit: concatenated

    script:
    """
    zcat ${read_files.join(' ')} | gzip > ${sample_id}_all_reads.fastq.gz
    """
}