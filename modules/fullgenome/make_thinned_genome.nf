process MAKE_THINNED_GENOME {
    tag "Create thinned genome for ${panel}"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0':
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0' }"

    storeDir "${params.genome_cache}/thinned_genomes/${panel}"

    input:
    val panel
    path genome
    path genome_fai
    path region_file

    output:
    path "thinned.fna", emit: fasta
    path "thinned.fna.fai", emit: fai
    path "thinned.fna.amb", emit: amb
    path "thinned.fna.ann", emit: ann
    path "thinned.fna.bwt", emit: bwt
    path "thinned.fna.pac", emit: pac
    path "thinned.fna.sa", emit: sa

    script:
    """
    samtools faidx ${genome} \$(cat ${region_file}) > thinned.fna
    bwa index thinned.fna
    samtools faidx thinned.fna
    """
}

process MAKE_THINNED_GENOME_MOUNT {
    tag "Create thinned genome for ${panel} (Mount)"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0':
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0' }"

    storeDir "${params.genome_cache}/thinned_genomes/${panel}"

    input:
    val panel
    val genome_path
    path region_file

    output:
    path "thinned.fna", emit: fasta
    path "thinned.fna.fai", emit: fai
    path "thinned.fna.amb", emit: amb
    path "thinned.fna.ann", emit: ann
    path "thinned.fna.bwt", emit: bwt
    path "thinned.fna.pac", emit: pac
    path "thinned.fna.sa", emit: sa

    script:
    """
    GENOME_FILE="${genome_path}/Otsh_v1.0/Otsh_v1.0.fna"

    # Verify file existence
    if [ ! -f "\$GENOME_FILE" ]; then
        echo "Error: Genome file not found at \$GENOME_FILE"
        echo "Please check 'full_genome_mount_path' parameter."
        exit 1
    fi
    
    # Run samtools on direct path
    samtools faidx \$GENOME_FILE \$(cat ${region_file}) > thinned.fna
    bwa index thinned.fna
    samtools faidx thinned.fna
    """
}
