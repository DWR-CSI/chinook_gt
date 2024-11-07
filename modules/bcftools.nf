process BCFTOOLS_MPILEUP {
    tag "MPILEUP on ${reference}"
    label 'process_high'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    
    publishDir "${params.outdir}/${params.project}/variants", mode: 'copy'

    input:
    tuple val(reference), path(bams), path(bais), path(ref)
    
    output:
    tuple val(reference), path("${params.project}_${reference}.vcf"), emit: vcf
    
    script:
    """
    bcftools mpileup -a AD,DP,INFO/AD \
        -B -q 20 -Q 20 \
        -f ${ref} \
        ${bams.join(' ')} \
        | bcftools call -v -m \
        | bcftools sort \
        | bcftools norm -m +any -f ${ref} \
        > ${params.project}_${reference}.vcf
    """
}