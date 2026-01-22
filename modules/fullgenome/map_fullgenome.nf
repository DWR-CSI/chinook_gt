process MAP_TO_FULL_GENOME {
    tag "Map ${sample_id} to full genome"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0':
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0' }"

    input:
    tuple val(sample_id), path(merged_reads)
    path genome
    path genome_index_files

    output:
    tuple val(sample_id), path("${sample_id}_fullg.bam"), path("${sample_id}_fullg.bam.bai"), emit: bam

    script:
    """
    # Write SAM to file first to avoid pipe buffering issues
    bwa mem \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}\\tLB:amplicon\\tPL:ILLUMINA\\tSM:${sample_id}" \
        ${params.full_genome_mount_path ? "${params.full_genome_mount_path}/Otsh_v1.0/Otsh_v1.0.fna" : genome} \
        ${merged_reads} \
        > ${sample_id}_fullg.sam

    # Convert to BAM and sort
    samtools view -u ${sample_id}_fullg.sam | \
        samtools sort -T ${sample_id}_tmp -O bam -o ${sample_id}_fullg.bam -

    samtools index ${sample_id}_fullg.bam
    """
}
