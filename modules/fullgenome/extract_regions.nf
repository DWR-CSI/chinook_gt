process EXTRACT_READS_FROM_REGIONS {
    tag "Extract ${sample_id} reads from target regions"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23--h96c455f_0':
        'quay.io/biocontainers/samtools:1.23--h96c455f_0' }"

    input:
    tuple val(sample_id), path(fullg_bam), path(fullg_bai)
    val panel
    path region_file

    output:
    tuple val(sample_id), val(panel), path("${sample_id}_${panel}_extracted.bam"), path("${sample_id}_${panel}_extracted.bam.bai"), emit: bam

    script:
    """
    # Extract reads from target regions to intermediate BAM
    samtools view -u -o ${sample_id}_${panel}_unsorted.bam ${fullg_bam} \$(cat ${region_file})

    # Sort the extracted reads
    samtools sort -T ${sample_id}_tmp -O bam -o ${sample_id}_${panel}_extracted.bam ${sample_id}_${panel}_unsorted.bam

    samtools index ${sample_id}_${panel}_extracted.bam
    """
}
