process TRIMMOMATIC {
    tag "Trimmomatic on ${sample_id}"
    label 'process_small'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    publishDir "${params.outdir}/${params.project}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), val(type), path(reads), path(adapter_file)
    val(trim_params)

    output:
    tuple val(sample_id), path("*_paired_*.fastq.gz"), emit: trimmed_paired
    tuple val(sample_id), path("*_unpaired_*.fastq.gz"), emit: trimmed_unpaired
    path "*.trimmomatic.log", emit: log

    script:
    def input_1 = reads[0]
    def input_2 = reads[1]
    def output_paired_1 = "${sample_id}_paired_R1.fastq.gz"
    def output_unpaired_1 = "${sample_id}_unpaired_R1.fastq.gz"
    def output_paired_2 = "${sample_id}_paired_R2.fastq.gz"
    def output_unpaired_2 = "${sample_id}_unpaired_R2.fastq.gz"
    def adapter_param = adapter_file ? "ILLUMINACLIP:${adapter_file}:2:30:10" : ''
    """
    trimmomatic PE \
        -threads ${task.cpus} \
        -phred33 \
        $input_1 $input_2 \
        $output_paired_1 $output_unpaired_1 \
        $output_paired_2 $output_unpaired_2 \
        $adapter_param \
        $trim_params \
        2> ${sample_id}.trimmomatic.log
    """
}

process TRIMMOMATIC_SINGLE {
    tag "Trimmomatic on ${sample_id}"
    label 'process_small'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    publishDir "${params.outdir}/${params.project}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), val(type), path(reads), path(adapter_file)
    val(trim_params)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: trimmed

    path "*.trimmomatic.log", emit: log

    script:
    def adapter_param = adapter_file ? "ILLUMINACLIP:${adapter_file}:2:30:10" : ''
    """
    trimmomatic SE \
        -threads ${task.cpus} \
        -phred33 \
        ${reads} \
        ${sample_id}_trimmed.fastq.gz \
        $adapter_param \
        $trim_params \
        2> ${sample_id}.trimmomatic.log
    """
}