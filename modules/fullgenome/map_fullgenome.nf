process MAP_TO_FULL_GENOME {
    tag "Map chunk of ${sample_ids instanceof List ? sample_ids.size() : 1} samples to full genome"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0':
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0' }"

    input:
    tuple val(sample_ids), path(merged_reads_list)
    path genome
    path genome_index_files

    output:
    tuple val(sample_ids), path("*_fullg.bam"), path("*_fullg.bam.bai"), emit: bam

    script:
    def samples = sample_ids instanceof List ? sample_ids : [sample_ids]
    def reads = merged_reads_list instanceof List ? merged_reads_list : [merged_reads_list]
    
    def bwa_cmds = ""
    for (int i = 0; i < samples.size(); i++) {
        def id = samples[i]
        def read = reads[i]
        bwa_cmds += """
        echo "Processing sample ${id}..."
        bwa mem \\
            -t ${task.cpus} \\
            -R "@RG\\\\tID:${id}\\\\tLB:amplicon\\\\tPL:ILLUMINA\\\\tSM:${id}" \\
            ${genome} \\
            ${read} \\
            > ${id}_fullg.sam
        
        samtools view -u ${id}_fullg.sam | \\
            samtools sort -T ${id}_tmp -O bam -o ${id}_fullg.bam -
        
        samtools index ${id}_fullg.bam
        rm ${id}_fullg.sam
        """
    }

    """
    ${bwa_cmds}
    """
}

process MAP_TO_FULL_GENOME_MOUNT {
    tag "Map chunk of ${sample_ids instanceof List ? sample_ids.size() : 1} samples to full genome (Mount)"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0':
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0' }"

    input:
    tuple val(sample_ids), path(merged_reads_list)
    val genome_path

    output:
    tuple val(sample_ids), path("*_fullg.bam"), path("*_fullg.bam.bai"), emit: bam

    script:
    def samples = sample_ids instanceof List ? sample_ids : [sample_ids]
    def reads = merged_reads_list instanceof List ? merged_reads_list : [merged_reads_list]

    def bwa_cmds = ""
    for (int i = 0; i < samples.size(); i++) {
        def id = samples[i]
        def read = reads[i]
        bwa_cmds += """
        echo "Processing sample ${id}..."
        if ! bwa mem \\
            -t ${task.cpus} \\
            -R "@RG\\\\tID:${id}\\\\tLB:amplicon\\\\tPL:ILLUMINA\\\\tSM:${id}" \\
            \$LOCAL_GENOME \\
            ${read} \\
            > ${id}_fullg.sam 2> bwa_${id}.stderr; then
            
            echo "BWA MEM FAILED for ${id}"
            echo "=== BWA STDERR BEGIN ==="
            cat bwa_${id}.stderr
            echo "=== BWA STDERR END ==="
            exit 1
        fi

        samtools view -u ${id}_fullg.sam | \\
            samtools sort -T ${id}_tmp -O bam -o ${id}_fullg.bam -
        
        samtools index ${id}_fullg.bam
        rm ${id}_fullg.sam
        rm bwa_${id}.stderr
        """
    }

    """
    GENOME_FILE="${genome_path}/Otsh_v1.0/Otsh_v1.0.fna"

    # Verify file existence
    if [ ! -f "\$GENOME_FILE" ]; then
        echo "Error: Genome file not found at \$GENOME_FILE"
        exit 1
    fi

    # Copy to local disk once
    echo "Copying genome to local disk..."
    cp "\$GENOME_FILE"* .
    LOCAL_GENOME="Otsh_v1.0.fna"

    ${bwa_cmds}
    """
}
