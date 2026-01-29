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
        ${genome} \
        ${merged_reads} \
        > ${sample_id}_fullg.sam

    # Convert to BAM and sort
    samtools view -u ${sample_id}_fullg.sam | \
        samtools sort -T ${sample_id}_tmp -O bam -o ${sample_id}_fullg.bam -

    samtools index ${sample_id}_fullg.bam
    """
}

process MAP_TO_FULL_GENOME_MOUNT {
    tag "Map ${sample_id} to full genome (Mount)"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0':
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0' }"

    input:
    tuple val(sample_id), path(merged_reads)
    val genome_path

    output:
    tuple val(sample_id), path("${sample_id}_fullg.bam"), path("${sample_id}_fullg.bam.bai"), emit: bam

    script:
    """
    GENOME_FILE="${genome_path}/Otsh_v1.0/Otsh_v1.0.fna"

    # Verify file existence to support clean failure
    if [ ! -f "\$GENOME_FILE" ]; then
        echo "Error: Genome file not found at \$GENOME_FILE"
        echo "Please check 'full_genome_mount_path' parameter."
        exit 1
    fi

    # Write SAM to file first to avoid pipe buffering issues
    # Capture stderr to identify mount/memory issues
    
    # Copy to local disk to avoid file share throttling
    echo "Copying genome to local disk..."
    cp "\$GENOME_FILE"* .
    LOCAL_GENOME="Otsh_v1.0.fna"

    if ! bwa mem \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}\\tLB:amplicon\\tPL:ILLUMINA\\tSM:${sample_id}" \
        \$LOCAL_GENOME \
        ${merged_reads} \
        > ${sample_id}_fullg.sam 2> bwa.stderr; then
        
        echo "BWA MEM FAILED"
        echo "=== BWA STDERR BEGIN ==="
        cat bwa.stderr
        echo "=== BWA STDERR END ==="
        exit 1
    fi

    # Convert to BAM and sort
    samtools view -u ${sample_id}_fullg.sam | \
        samtools sort -T ${sample_id}_tmp -O bam -o ${sample_id}_fullg.bam -

    samtools index ${sample_id}_fullg.bam
    """
}
