process BAM_TO_FASTQ {
    tag "Convert ${sample_id} extracted BAM to FASTQ"
    label 'process_small'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3':
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    tuple val(sample_id), val(panel), path(extracted_bam), path(extracted_bai)

    output:
    tuple val(sample_id), val(panel), path("${sample_id}_${panel}_extracted.fastq.gz"), emit: fastq

    script:
    """
    bedtools bamtofastq -i ${extracted_bam} -fq ${sample_id}_${panel}_extracted.fastq
    gzip ${sample_id}_${panel}_extracted.fastq
    """
}
