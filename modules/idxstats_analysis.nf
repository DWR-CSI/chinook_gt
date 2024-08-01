process ANALYZE_IDXSTATS {
    tag "Analyzing idxstats"
    label 'process_small'
    publishDir "${params.outdir}/${params.project}/idxstats_analysis", mode: 'copy'
    container "docker.io/bnguyen29/r-rubias:1.0.2"

    input:
    path idxstats_files

    output:
    path "*.{pdf,jpg,png}", emit: plots
    path "reads_matrix.txt", emit: matrix
    
    script:
    """
    Rscript $projectDir/bin/analyze_idxstats.R . $PWD ${params.n_loci}
    """
}