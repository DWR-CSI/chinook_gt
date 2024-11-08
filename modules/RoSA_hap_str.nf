process GREB_HAPSTR {
    tag 'Extracting RoSA info via hapstr'
    container 'quay.io/biocontainers/mulled-v2-9a0dd936806dc33b863c3b3f851665f40f2af214:5aaef6d70fcd32819000b9ac3c8ce065fbcaac8b-0'
    label 'process_medium'

    publishDir "${params.outdir}/${params.project}/rosa", mode: 'copy'

    input:
    tuple(val(reference),path(filtered_vcf))
    path(allele_key)

    output:
    path "*_ots28_report.tsv", emit: ots28_report
    path("${params.project}_RoSA_hapstrs.tsv"), emit: rosa_hapstrs

    script:
    """
    rosa_hapstr.R ${params.project} $filtered_vcf $allele_key
    """
}