process BCFTOOLS_MPILEUP {
    tag "MPILEUP on ${reference}"
    label 'process_high'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    
    publishDir "${params.outdir}/${params.project}/variants", mode: 'copy'

    input:
    tuple val(reference), path(bams), path(bais), path(ref)
    
    output:
    tuple val(reference), path("${params.project}_${reference}.vcf"), emit: vcf
    tuple val(reference), path("${params.project}_${reference}.filtered.vcf"), emit: filtered_vcf, optional: true
    
    script:
    def filter_cmd = (reference == "transition" || reference == "full") ? 
        """
        bcftools view -i 'QUAL>=20 && FORMAT/DP>=5' --exclude-types indels \
            ${params.project}_${reference}.vcf \
            -O v \
            -o ${params.project}_${reference}.filtered.vcf
        """ : 
        ""

    """
    bcftools mpileup -a AD,DP,INFO/AD \
        -B -q 20 -Q 20 \
        -f ${ref} \
        ${bams.join(' ')} \
        | bcftools call -v -m \
        | bcftools sort \
        | bcftools norm -m +any -f ${ref} \
        > ${params.project}_${reference}.vcf
    
    ${filter_cmd}
    """
}