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
    set -euo pipefail

    # write header from first file
    first=\$(ls *.counts.tsv | head -n 1)
    head -n 1 \$first > combined_dimer_counts.tsv

    # append all files
    for f in *.counts.tsv; do
        tail -n +2 \$f >> combined_dimer_counts.tsv
    done
    """
}
