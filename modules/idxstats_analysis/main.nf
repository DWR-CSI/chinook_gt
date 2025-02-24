process ANALYZE_IDXSTATS {
    tag "Analyzing idxstats and loci distribution"
    label 'process_small'
    publishDir "${params.outdir}/${params.project}/idxstats_analysis", mode: 'copy'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-70a811298b15e02937b3477ebadce5062ab47a53:8e6675e8830568a8642779e4ba23787ea188a067-0':
        'quay.io/biocontainers/mulled-v2-70a811298b15e02937b3477ebadce5062ab47a53:8e6675e8830568a8642779e4ba23787ea188a067-0' }"

    input:
    path idxstats_files

    output:
    path "idxstats_plots/*.{pdf,jpg,png}", emit: idxstats_plots
    path "loci_qc/*.{pdf,csv,txt}", emit: loci_qc
    path "reads_matrix.txt", emit: matrix
    
    script:
    """
    # Create output directories
    mkdir -p idxstats_plots loci_qc
    
    # Run original idxstats analysis
    analyze_idxstats.R . idxstats_plots ${params.n_loci}
    
    # Run additional loci QC analysis
    check_loci.R idxstats_plots/reads_matrix.txt loci_qc
    """
}