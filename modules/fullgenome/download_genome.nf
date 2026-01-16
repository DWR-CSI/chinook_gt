process DOWNLOAD_AND_INDEX_GENOME {
    tag "Download and index Otsh_v1.0 genome"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-30758c5e2be407d6bbf2570c832fb245e8488634:33e45705f656cfb3e5268b3d075c838c719a0fa2':
        'quay.io/biocontainers/mulled-v2-30758c5e2be407d6bbf2570c832fb245e8488634:33e45705f656cfb3e5268b3d075c838c719a0fa2' }"

    storeDir "${params.genome_cache}/Otsh_v1.0"

    output:
    path "Otsh_v1.0.fna", emit: genome
    path "Otsh_v1.0.fna.fai", emit: fai
    path "Otsh_v1.0.fna.amb", emit: amb
    path "Otsh_v1.0.fna.ann", emit: ann
    path "Otsh_v1.0.fna.bwt", emit: bwt
    path "Otsh_v1.0.fna.pac", emit: pac
    path "Otsh_v1.0.fna.sa", emit: sa

    script:
    """
    wget -O Otsh_v1.0.fna.gz "${params.genome_url}"
    gunzip Otsh_v1.0.fna.gz
    bwa index Otsh_v1.0.fna
    samtools faidx Otsh_v1.0.fna
    """
}
