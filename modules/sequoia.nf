process RUN_SEQUOIA {
    tag "Parentage Analysis: sequoia"
    // container "docker.io/bnguyen29/r-rubias:1.0.4" // TO-DO: needs to be replaced with sequoia container when available
    label 'process_high'
    publishDir "${params.outdir}/${params.project}/sequoia", mode: 'copy'

    input:
    path(parent_geno)
    path(parent_lifehistory)
    path(offspring_geno)
    val(offspring_BY)
    val(offspring_minBY)
    val(offspring_maxBY)
    val(offspring_max_age)

    output:
    path "*_sequoia_output.rds", emit: rds_output
    path "*_parentage_results.txt", emit: parentage
    path "*_sequoia_genotype_matrix.txt", emit: genotype_matrix
    path "*_sequoia_lifehistory_data.txt", emit: lh
    path "*_allele_dictionary.txt", emit: allele_dict

    script:
    """
    run_sequoia.R ${params.sequoia_mode} ${parent_geno} ${parent_lifehistory} ${offspring_geno} ${offspring_BY} ${offspring_minBY} ${offspring_maxBY} ${offspring_max_age} ${params.project} ${params.sequoia_missing_threshold}
    """
}
