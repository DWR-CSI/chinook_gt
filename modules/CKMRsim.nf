process CKMR_GENO_LONG {
    tag "CKMR Genotype Long Format Conversion"
    container "docker.io/bnguyen29/ckmrsim:latest"
    label 'process_medium'
    publishDir "${params.outdir}/${params.project}/CKMRsim", mode: 'copy'

    input:
    path(geno_wide)
    path(parents_geno_wide)
    path(allele_frequencies)
    val(logl_threshold)

    output:
    path "*_PO_results.tsv", emit: PO_results
    path "*_parent_genotypes_long.tsv", emit: parent_genotypes_long
    path "*_offspring_genotypes_long.tsv", emit: offspring_genotypes_long

    script:
    """
    export LOCI_REMOVAL_REGEX='${params.loci_to_remove}'
    run_CKMR.R ${geno_wide} ${parents_geno_wide} ${logl_threshold} ${params.project} ${params.CKMR_min_loci} ${allele_frequencies}
    """
}