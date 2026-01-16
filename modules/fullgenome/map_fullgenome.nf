process MAP_TO_FULL_GENOME {
    tag "Map ${sample_id} to full genome"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-30758c5e2be407d6bbf2570c832fb245e8488634:33e45705f656cfb3e5268b3d075c838c719a0fa2':
        'quay.io/biocontainers/mulled-v2-30758c5e2be407d6bbf2570c832fb245e8488634:33e45705f656cfb3e5268b3d075c838c719a0fa2' }"

    input:
    tuple val(sample_id), path(merged_reads)
    path genome
    path genome_index_files

    output:
    tuple val(sample_id), path("${sample_id}_fullg.bam"), path("${sample_id}_fullg.bam.bai"), emit: bam

    script:
    """
    bwa mem \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}\\tLB:amplicon\\tPL:ILLUMINA\\tSM:${sample_id}" \
        ${genome} \
        ${merged_reads} \
    | samtools view -u - \
    | samtools sort -T ${sample_id}_tmp -O bam -o ${sample_id}_fullg.bam -

    samtools index ${sample_id}_fullg.bam
    """
}
