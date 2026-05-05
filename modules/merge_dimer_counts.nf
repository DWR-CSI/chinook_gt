process MERGE_DIMER_COUNTS {
    tag "Merge dimer counts"
    label 'process_single'

    publishDir "${params.outdir}/${params.project}/dimer_counts", mode: 'copy'

    input:
    path(count_files)

    output:
    path "combined_dimer_counts.tsv"

    script:
    """
    # write header once
    head -n 1 \$(ls ${count_files[0]}) > combined_dimer_counts.tsv

    # append all rows (skip headers)
    for f in ${count_files}; do
        tail -n +2 \$f >> combined_dimer_counts.tsv
    done
    """
}
