process DIMER_COUNTS {
    tag "Combine FLASH2 dimer counts"
    label 'process_single'
    publishDir "${params.outdir}/${params.project}/dimer_counts", mode: 'copy'

    input:
    path(count_files)

    output:
    path "combined_dimer_counts.tsv", emit: combined

    script:
    """
    # Write header
    echo -e "sample_id\\ttotal_reads\\tdimer_reads\\tpercent_dimer" > combined_dimer_counts.tsv
    
    # Concatenate all count files, skipping headers
    for file in ${count_files}; do
        tail -n +2 "\$file" >> combined_dimer_counts.tsv
    done
    """
}
