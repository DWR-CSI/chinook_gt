process BWA_MEM {
    tag "Mapping ${sample_id} to $ref_name"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bwa:0.7.18--he4a0461_1':
    'quay.io/biocontainers/bwa:0.7.18--he4a0461_1' }"

    publishDir "${params.outdir}/${params.project}/bwa_mem/${ref_name}", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_reads), val(ref_name), path(ref_and_index_files), path(vcf)
    // VCF file is not used right now.

    output:
    tuple val(sample_id), val(ref_name), path("${sample_id}_${ref_name}_aln.sam"), emit: aligned_sam

    script:
    def reference_fasta = ref_and_index_files.find { it.name =~ /\.(fasta|fna)$/ }
    """
    bwa mem \
        -M \
        -v 3 \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}\\tLB:amplicon\\tPL:ILLUMINA\\tSM:${sample_id}" \
        "${reference_fasta}" \
        "${merged_reads}" \
        > "${sample_id}_${ref_name}_aln.sam"
    """
}