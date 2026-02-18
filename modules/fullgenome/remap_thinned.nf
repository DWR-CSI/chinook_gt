process REMAP_TO_THINNED_GENOME {
    tag "Remap ${sample_id} to thinned ${panel} genome"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bwa:0.7.19--h577a1d6_1':
    'quay.io/biocontainers/bwa:0.7.19--h577a1d6_1' }"

    input:
    tuple val(sample_id), val(panel), path(extracted_fastq)
    path thinned_fasta
    path thinned_index_files

    output:
    tuple val(sample_id), val(panel), path("${sample_id}_${panel}_aln.sam"), emit: aligned_sam

    script:
    """
    bwa mem \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}\\tLB:amplicon\\tPL:ILLUMINA\\tSM:${sample_id}" \
        ${thinned_fasta} \
        ${extracted_fastq} \
        > ${sample_id}_${panel}_aln.sam
    """
}
