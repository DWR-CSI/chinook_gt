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
    samtools faidx ${params.full_genome_mount_path ? "${params.full_genome_mount_path}/reference/genomes/Otsh_v1.0/Otsh_v1.0.fna" : genome} \$(cat ${region_file}) > thinned.fna
    bwa index thinned.fna
    samtools faidx thinned.fna
    """
}
