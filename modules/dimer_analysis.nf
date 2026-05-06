process DIMER_ANALYSIS {
    tag "Dimer analysis on ${sample_id}"
    label 'process_small'
    container 'quay.io/biocontainers/flash2:2.2.00--ha92aebf_1'

    input:
    tuple val(sample_id), val(read_type), path(reads)
    val min_overlap
    val min_outie_overlap
    val max_overlap
    
    output:
    path "${sample_id}.counts.tsv", emit: counts
    
    script:
    def input_1 = reads[0]
    def input_2 = reads[1]

    """
    flash2 \
        -m ${min_overlap} \
        --min-overlap-outie ${min_outie_overlap} \
        -M ${max_overlap} \
        -O \
        -z \
        -o ${sample_id} \
        -d . \
        ${input_1} ${input_2} \
        > ${sample_id}.flash2.log 2>&1

	//compute metrics and write ONE LINE ONLY (no header)
    zcat ${sample_id}.extendedFrags.fastq.gz | \
    awk -v sid="${sample_id}" '
        NR % 4 == 2 {
            total++
            if (length(\$0) <= 50) short++
        }
        END {
            pct = (total > 0) ? (short / total) * 100 : 0
            printf "%s\\t%d\\t%d\\t%.2f\\n", sid, total, short, pct
        }
    ' > ${sample_id}.counts.tsv

	rm -f ${sample_id}.extendedFrags.fastq.gz \
      	${sample_id}.notCombined_*.fastq.gz

    """
}
