process CKMR_PO {
    tag "CKMR Parent-Offspring Analysis"
    container "docker.io/bnguyen29/ckmrsim:latest"
    label 'process_medium'
    publishDir "${params.outdir}/${params.project}/CKMRsim", mode: 'copy'

    input:
    path(geno_wide)
    path(parents_geno_wide)
    path(extra_genos_allele_freqs)

    output:
    path "*_PO_results.tsv", emit: PO_results
    path "*_PO_ckmr.rds", emit: PO_ckmr
    path "*_offspring_genotypes_long.rds", emit: offspring_genotypes_long
    path "*_parent_genotypes_long.rds", emit: parent_genotypes_long

    script:
    """
    export LOCI_REMOVAL_REGEX='${params.loci_to_remove}'
    run_CKMR.R ${geno_wide} ${parents_geno_wide} ${params.CKMR_logl_threshold} ${params.project} ${params.CKMR_min_loci} ${extra_genos_allele_freqs}
    """
}