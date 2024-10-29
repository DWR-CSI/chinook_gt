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

def resolveReferences(params) {
    def reference_files = []
    
    if (params.reference) {
        // User provided explicit reference files
        if (params.reference instanceof List) {
            reference_files = params.reference
        } else {
            reference_files = [params.reference]
        }
        log.info "Using user-specified reference file(s): ${reference_files}"
    } else if (params.panel) {
        // Handle panel-based reference selection
        switch(params.panel.toLowerCase()) {
            case 'transition':
                reference_files = [
                    "$projectDir/data/targets/transition/transition.fasta",
                    //"$projectDir/data/targets/LFAR/LFAR.fna", // Commented out. Not used in DWR primers?
                    "$projectDir/data/targets/WRAP/WRAP.fna"
                ]
                log.info "Using transition panel references"
                break
            case 'full':
                reference_files = [
                    "$projectDir/data/targets/full/full.fna",
                    "$projectDir/data/targets/full_VGLL3Six6LFARWRAP/VGLL3Six6LFARWRAP.fna"
                ]
                log.info "Using full panel reference"
                break
            default:
                error "Unrecognized panel type: ${params.panel}. Supported values are 'transition' or 'full'"
        }
        
        // Verify all reference files exist
        reference_files.each { ref ->
            if (!file(ref).exists()) {
                error "Reference file not found: ${ref}"
            }
        }
        
        // Verify corresponding VCF files exist
        reference_files.each { ref ->
            def basename = file(ref).simpleName
            def vcf = file("${projectDir}/data/VCFs/${basename}.vcf")
            if (!vcf.exists()) {
                error "VCF file not found for reference: ${basename}. Expected: ${vcf}"
            }
        }
    } else {
        error "Neither reference nor panel type specified. Please provide either --reference or --panel (transition/full)"
    }
    
    return reference_files
}
    
workflow {
    log.info """
    ==============================================
    DWR-CSI/chinook_gt Pipeline
    ==============================================
    Panel Type : ${params.panel ?: 'Not specified'}
    References : ${params.reference ? 'User specified' : 'Auto-selected'}
    """
    // Resolve and validate references
    reference_files = resolveReferences(params)
    // Create reference channel with associated files
    Channel
        .fromList(reference_files)
        .flatMap { ref ->
            // For each reference, check all required index files
            def ref_file = file(ref)
            def required_extensions = ['', '.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']
            def missing_files = []
            
            def ref_files = required_extensions.collect { ext ->
                def f = file("${ref}${ext}")
                if (!f.exists()) {
                    missing_files << "${ref}${ext}"
                }
                return f
            }
            
            if (missing_files) {
                log.warn "Missing index files for ${ref}: ${missing_files}"
                log.info "Running BWA index to generate missing files..."
                // You might want to add an INDEX_REFERENCE process here
            }
            
            return ref_files
        }
        .collect()
        .map { files -> 
            files.groupBy { it.simpleName }
        }
        .flatMap { refMap -> 
            refMap.collect { refName, refFiles ->
                // Also locate the corresponding VCF file
                def vcf = file("${projectDir}/data/VCFs/${refName}.vcf")
                tuple(refName, refFiles, vcf)
            }
        }
        .set { reference_ch }
    
    // Define input channels
    Channel
        .fromFilePairs(params.input, checkIfExists: true)
        .set { ch_input_fastq_pairs }
        
    ch_input_fastq = Channel
        .fromPath(params.input, checkIfExists: true)
        .collect()
    
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

    BWA_MEM.out.aligned_sam
        .groupTuple(by: 1)
        .map { sample_id, ref_name, sam_files ->
            // Search for VCF file in data/VCFs directory of projectDir
            // VCF file should have the same basename as the reference files.
            def vcfFile = file("${projectDir}/data/VCFs/${ref_name}.vcf")
            if (!vcfFile.exists()) {
                error "VCF file not found for reference: ${ref_name}. Place a file named ${ref_name}.vcf in the data/VCFs of project directory."
            }
            tuple(sample_id, ref_name, vcfFile, sam_files)
        }
        .set { grouped_sam_files }
    //grouped_sam_files.view { println "Debug: Grouped SAM files: $it" }

    // Now you can use grouped_sam_files in your next process
    GEN_MHP_SAMPLE_SHEET(grouped_sam_files)

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
    PREP_MHP_RDS(GEN_MHP_SAMPLE_SHEET.out.mhp_samplesheet)
    GEN_HAPS(PREP_MHP_RDS.out.rds)
    panel_branched_haps = GEN_HAPS.out.haps
        .branch { 
            main: it[0] =~ /transition|full/
            other: true
            }
    HAP2GENO(panel_branched_haps.main, locus_index_ch)
    CHECK_FILE_UPDATE(HAP2GENO.out.new_index, locus_index_ch) | view

    // Run Structure and Rubias analyses
    STRUC_PARAMS(ots28_baseline_ch, HAP2GENO.out.numgeno_OTS28)
    STRUCTURE(STRUC_PARAMS.out.structure_input, STRUC_PARAMS.out.m_params, STRUC_PARAMS.out.e_params)
    STRUCTURE_ROSA_REPORT(STRUCTURE.out.structure_output, STRUC_PARAMS.out.structure_input, params.ots28_missing_threshold)
    RUN_RUBIAS(STRUCTURE_ROSA_REPORT.out.ots28_report, HAP2GENO.out.numgeno, baseline_ch)
}
