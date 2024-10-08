process SAMTOOLS {
    tag "Create BAM and IDXSTATS: ${sample_id}"
    label 'process_small'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    publishDir "${params.outdir}/${params.project}/samtools/${reference}", mode: 'copy', pattern: '*_sorted.bam*'
    publishDir "${params.outdir}/${params.project}/samtools/idxstats/${reference}", mode: 'copy', pattern: '*_idxstats.txt'

    input:
    tuple val(sample_id), val(reference), path(aligned_sam)

    output:
    tuple val(sample_id), val(reference), path("${sample_id}_sorted.bam"), emit: sorted_bam
    tuple val(sample_id), val(reference), path("${sample_id}_sorted.bam.bai"), emit: sorted_bam_index
    tuple val(sample_id), val(reference), path("${sample_id}_idxstats.txt"), emit: idxstats

    script:
    """
    samtools view -bS "${aligned_sam}" | samtools sort -o "${sample_id}_sorted.bam" -
    samtools index "${sample_id}_sorted.bam"
    samtools idxstats "${sample_id}_sorted.bam" > "${sample_id}_idxstats.txt"
    """

}