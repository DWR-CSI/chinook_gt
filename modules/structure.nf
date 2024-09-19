process STRUC_PARAMS {
    tag "Prepare STRUCTURE param files"
    container 'ubuntu:latest'
    label 'process_small'

    publishDir "${params.outdir}/${params.project}/structure", mode: 'copy'

    input:
    path ots28_baseline
    path ots28_num_genotypes

    output:
    path "mainparams.txt", emit: m_params
    path "extraparams.txt", emit: e_params
    path "structure_data.txt", emit: structure_input

    shell:
    """
    #!/bin/bash
    cat !{ots28_baseline} > structure_data.txt
    tail -n +2 !{ots28_num_genotypes} >> structure_data.txt
    N_LINES=\$(wc -l < structure_data.txt)
    N_INDIVS=\$((N_LINES - 1))

    cat <<EOF > mainparams.txt
    #define OUTFILE !{params.project}_RoSA_structure.txt
    #define INFILE structure_data.txt
    #define NUMINDS \$N_INDIVS
    #define NUMLOCI 7
    #define LABEL 1 
    #define POPDATA 1 
    #define POPFLAG 1 
    #define LOCDATA 0 
    #define PHENOTYPE 0 
    #define MARKERNAMES 1 
    #define MAPDISTANCES 0 
    #define ONEROWPERIND 1 
    #define PHASEINFO 0 
    #define PHASED 0 
    #define RECESSIVEALLELES 0 
    #define EXTRACOLS 0
    #define MISSING -9
    #define PLOIDY 2
    #define MAXPOPS 2
    #define BURNIN 25000
    #define NUMREPS 250000


    #define NOADMIX 0
    #define LINKAGE 0
    #define USEPOPINFO 0

    #define LOCPRIOR 0
    #define INFERALPHA 1
    #define ALPHA 1.0
    #define POPALPHAS 0 
    #define UNIFPRIORALPHA 1 
    #define ALPHAMAX 10.0
    #define ALPHAPROPSD 0.025


    #define FREQSCORR 1 
    #define ONEFST 0
    #define FPRIORMEAN 0.01
    #define FPRIORSD 0.05


    #define INFERLAMBDA 0 
    #define LAMBDA 1.0
    #define COMPUTEPROB 1 
    #define PFROMPOPFLAGONLY 0 
    #define ANCESTDIST 0 
    #define STARTATPOPINFO 0 
    #define METROFREQ 10


    #define UPDATEFREQ 1
    EOF

    touch extraparams.txt
    """
}

process STRUCTURE {
    tag "Run STRUCTURE"
    label 'process_small'
    container 'docker.io/rtibiocloud/structure:v2.3.4_f2d7e82'
    publishDir "${params.outdir}/${params.project}/structure", mode: 'copy'

    input:
    path genotypes_input
    path m_params
    path e_params
    

    output:
    path "*_RoSA_structure*", emit: structure_output

    script:
    """
    structure -m $m_params -e $e_params
    """
}