process MERGE_DIMER_COUNTS {
    input:
    path(count_files)

    output:
    path "combined_dimer_counts.tsv"

    script:
    """
    set -euo pipefail

    first=\$(ls *.counts.tsv | head -n 1)
    head -n 1 \$first > combined_dimer_counts.tsv

    for f in *.counts.tsv; do
        tail -n +2 \$f >> combined_dimer_counts.tsv
    done
    """
}
