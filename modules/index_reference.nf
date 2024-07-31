process INDEX_REFERENCE {
    tag "Indexing reference genome"
    label 'process_high'
    container 'docker.io/biocontainers/bwa:v0.7.17_cv1'
    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*")

    script:
    """
    bwa index ${reference}
    """
}