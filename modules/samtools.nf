process SAMTOOLS {
    tag "SAMTOOLS on ${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    publishDir "${params.outdir}/${params.project}/samtools", mode: 'copy'

    input:
    tuple val(sample_id), path(aligned_sam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bam
    tuple val(sample_id), path("${sample_id}_sorted.bam.bai"), emit: sorted_bam_index
    tuple val(sample_id), path("${sample_id}_idxstats.txt"), emit: idxstats
    path "${sample_id}_samtools.log", emit: log

    script:
    """
    samtools view -bS "${aligned_sam}" | samtools sort -o "${sample_id}_sorted.bam" -
    samtools index "${sample_id}_sorted.bam"
    samtools idxstats "${sample_id}_sorted.bam" > "${sample_id}_idxstats.txt"
    """

}