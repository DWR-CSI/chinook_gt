process CONCAT_READS {
    tag "Concatenating reads for ${sample_id}"
    container 'docker.io/nfcore/base:2.1'
    label 'process_small'
    
    publishDir "${params.outdir}/${params.project}/concatenated", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    tuple val(sample_id), path("${sample_id}_all_reads.fastq.gz"), emit: concatenated

    script:
    """
    # Handle potential file path issues on Azure by using explicit file listing
    for file in ${read_files.join(' ')}; do
        if [[ -f "\$file" ]]; then
            cat "\$file" >> ${sample_id}_all_reads.fastq.gz
        else
            echo "Warning: File \$file not found" >&2
        fi
    done
    
    # Ensure output file exists even if no input files were found
    if [[ ! -f ${sample_id}_all_reads.fastq.gz ]]; then
        touch ${sample_id}_all_reads.fastq.gz
    fi
    """
}