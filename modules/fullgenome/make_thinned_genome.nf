process MAKE_THINNED_GENOME {
    tag "Create thinned genome for ${panel}"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-30758c5e2be407d6bbf2570c832fb245e8488634:33e45705f656cfb3e5268b3d075c838c719a0fa2':
        'quay.io/biocontainers/mulled-v2-30758c5e2be407d6bbf2570c832fb245e8488634:33e45705f656cfb3e5268b3d075c838c719a0fa2' }"

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
