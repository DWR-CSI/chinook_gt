process BWA_MEM {
    tag "BWA MEM on ${sample_id}"
    label 'process_high'
    container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'

    publishDir "${params.outdir}/${params.project}/bwa_mem", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_reads)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_aln.sam"), emit: aligned_sam
    path "${sample_id}_bwa.log", emit: log

    script:
    """
    bwa mem \
        -M \
        -v 3 \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}\\tLB:amplicon\\tPL:ILLUMINA\\tSM:${sample_id}" \
        "${reference}" \
        "${merged_reads}" \
        > "${sample_id}_aln.sam" \
        2> "${sample_id}_bwa.log"
    """
}