#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DWR-CSI/chinook_gt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/DWR-CSI/chinook_gt
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Import modules
include { FASTQC } from './modules/fastqc'
include { TRIMMOMATIC } from './modules/trimmomatic'
include { FLASH2 } from './modules/flash2'
include { BWA_MEM } from './modules/bwa_mem'
include { SAMTOOLS } from './modules/samtools'
include { MULTIQC } from './modules/multiqc'
include { INDEX_REFERENCE } from './modules/index_reference'
include { ANALYZE_IDXSTATS } from './modules/idxstats_analysis'
include { GEN_MHP_SAMPLE_SHEET; PREP_MHP_RDS; GEN_HAPS; HAP2GENO; CHECK_FILE_UPDATE } from './modules/microhaplot.nf'
include { RUN_RUBIAS } from './modules/rubias.nf'
include { STRUC_PARAMS; STRUCTURE } from './modules/structure.nf'
include { STRUCTURE_ROSA_REPORT } from './modules/rosa.nf'

// Functions


    
workflow {
    if (params.reference == null) {
        if (params.panel == null) {
            error "No reference or panel file provided. Exiting."
        } else {
            if (params.panel.toLowerCase() == "transition") { // Transition panel defaults
                params.reference = "$projectDir/data/transition_target/transition.fasta"
                if ((params.addtl_target == null) && params.addtl_vcf == null) {
                    params.addtl_target = ["$projectDir/data/targets/LFAR/LFAR.fna", "$projectDir/data/targets/WRAP/WRAP.fna"]
                    params.addtl_vcf = ["$projectDir/data/targets/LFAR/LFAR.vcf", "$projectDir/data/targets/WRAP/WRAP.vcf"]
                }
            }
        }
    } else {
        if (params.panel != null) {
            "Both reference and panel specified. Looking for reference file."
        }
    }
    // Define input channels
    Channel // paired
        .fromFilePairs(params.input, checkIfExists: true)
        .set { ch_input_fastq_pairs }
    ch_input_fastq = Channel // all fastq files for FastQC
        .fromPath(params.input, checkIfExists: true)
        .collect()
    Channel // Load main reference file
        .fromPath("${params.reference}{,.amb,.ann,.bwt,.fai,.pac,.sa}")
                       .collect()
                       .map { files -> tuple(files[0].simpleName, files) }
        .set { reference_ch }
    view(reference_ch)
    
    adapters_ch = channel.fromPath(params.adapter_file)
    locus_index_ch = channel.fromPath(params.locus_index)
    baseline_ch = channel.fromPath(params.baseline)
    ots28_baseline_ch = channel.fromPath(params.ots28_baseline)
    FASTQC(ch_input_fastq)
    ch_fastq_adapters_combined = ch_input_fastq_pairs.combine(adapters_ch)
    TRIMMOMATIC(ch_fastq_adapters_combined, params.trim_params)
    FLASH2(TRIMMOMATIC.out.trimmed_paired, params.min_overlap, params.max_overlap)
    bwa_input = FLASH2.out.merged.combine(reference_ch)
    BWA_MEM(bwa_input)
    SAMTOOLS(BWA_MEM.out.aligned_sam)
    // Collect all idxstats files
    idxstats_files_sorted = SAMTOOLS.out.idxstats
        .branch { 
            main: it[1] =~ /transition|full/
            other: true
            }

    idxstats_main = idxstats_files_sorted.main
        .map { it -> it[2] }  // Extract just the file path from the tuple
        .collect()            // Collect all files into a single list

    // Run the analysis
    ANALYZE_IDXSTATS(idxstats_main)

    branched_sams = BWA_MEM.out.aligned_sam
        .branch { 
            main: it[1] =~ /transition|full/
            other: true
            }

    GEN_MHP_SAMPLE_SHEET(branched_sams.main.map { it -> it[2] }.collect())

    // Run additional marker sets
    /*
    if (params.addtl_target && params.addtl_vcf) {
        addtl_vcf = Channel.fromPath(params.addtl_vcf, checkIfExists: true)
        Channel
            .fromPath(params.addtl_target, checkIfExists: true)
            .merge(addtl_vcf)
            .set { addtl_target_vcf }
        bwa_addtl_input = FLASH2.out.merged.combine(addtl_target_vcf)
    } else {
        addtl_target_vcf = Channel.empty()
    }
    */
    

    
    // Collect all QC files
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.fastqc_results)
    ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log)
    ch_multiqc_files = ch_multiqc_files.mix(FLASH2.out.log)
    //ch_multiqc_files = ch_multiqc_files.mix(BWA_MEM.out.log)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS.out.idxstats.collect{it[2]})
    //ch_multiqc_files.view { println "Debug: MultiQC input file: $it" }
    MULTIQC(ch_multiqc_files.collect())

    // Run the microhaplotype analysis
    vcf_file_ch = Channel.fromPath(params.baseline_vcf)
    PREP_MHP_RDS(GEN_MHP_SAMPLE_SHEET.out.mhp_samplesheet, vcf_file_ch, branched_sams.main.map { it -> it[2] }.collect())
    GEN_HAPS(PREP_MHP_RDS.out.rds)
    HAP2GENO(GEN_HAPS.out.haps, locus_index_ch)
    CHECK_FILE_UPDATE(HAP2GENO.out.new_index, locus_index_ch) | view

    // Run Structure and Rubias analyses
    STRUC_PARAMS(ots28_baseline_ch, HAP2GENO.out.numgeno_OTS28)
    STRUCTURE(STRUC_PARAMS.out.structure_input, STRUC_PARAMS.out.m_params, STRUC_PARAMS.out.e_params)
    STRUCTURE_ROSA_REPORT(STRUCTURE.out.structure_output, STRUC_PARAMS.out.structure_input, params.ots28_missing_threshold)
    RUN_RUBIAS(HAP2GENO.out.numgeno, baseline_ch)
}

