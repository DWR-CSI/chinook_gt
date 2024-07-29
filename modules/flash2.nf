process FLASH2 {
    tag "FLASH2 on ${sample_id}"
    label 'process_small'
    container 'quay.io/biocontainers/flash2:2.2.00--ha92aebf_1'

    publishDir "${params.outdir}/${params.project}/flash2", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val min_overlap
    val max_overlap
    
    output:
    tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz"), emit: merged
    tuple val(sample_id), path("${sample_id}.notCombined_*.fastq.gz"), emit: unmerged
    path "${sample_id}.flash2.log", emit: log
    
    script:
    def input_1 = reads[0]
    def input_2 = reads[1]
    """
    flash2 \
        -m ${min_overlap} \
        -M ${max_overlap} \
        -O \
        -z \
        -o ${sample_id} \
        -d . \
        ${input_1} ${input_2} \
        > ${sample_id}.flash2.log 2>&1
    """
}