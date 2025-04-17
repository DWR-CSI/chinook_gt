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

// Set default values for new rubias parameters
params.ots28_missing_threshold = params.ots28_missing_threshold ?: 0.5
params.gsi_missing_threshold = params.gsi_missing_threshold ?: 0.6
params.pofz_threshold = params.pofz_threshold ?: 0.8


// Import modules
include { FASTQC } from './modules/fastqc'
include { TRIMMOMATIC; TRIMMOMATIC_SINGLE } from './modules/trimmomatic'
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
include { BCFTOOLS_MPILEUP } from './modules/bcftools.nf'
include { GREB_HAPSTR } from './modules/RoSA_hap_str.nf'

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
                    //"$projectDir/data/targets/LFAR/LFAR.fasta", // Commented out. Not used in DWR primers?
                    "$projectDir/data/targets/WRAP/WRAP.fasta"
                ]
                log.info "Using transition panel references"
                break
            case 'full':
                reference_files = [
                    "$projectDir/data/targets/full/full.fna"//,
                    //"$projectDir/data/targets/full_VGLL3Six6LFARWRAP/VGLL3Six6LFARWRAP.fna"
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
    
    Rubias parameters:
    OTS28 Missing Threshold : ${params.ots28_missing_threshold}
    GSI Missing Threshold   : ${params.gsi_missing_threshold}
    PofZ Threshold          : ${params.pofz_threshold}
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
                log.info "Please run BWA index to generate missing files..."
                error "Fatal error: Missing required index files for ${ref}. Please ensure all necessary index files are present."
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
    adapters_ch = channel.fromPath(params.adapter_file)

    def locus_index_file
    if (params.containsKey('locus_index')) {
        locus_index_file = params.locus_index
    } else {
        locus_index_file = "${projectDir}/data/indices/transition_loci_indices_20241104.csv"
        log.info "No locus index provided, using default: ${locus_index_file}"
    }
    locus_index_ch = channel.fromPath(locus_index_file, checkIfExists: true)

    def baseline_file
    if (params.containsKey('baseline')) {
        baseline_file = params.baseline
    } else {
        baseline_file = "${projectDir}/data/baselines/full/SWFSC-chinook-reference-baseline-CV.csv"
        log.info "No baseline provided, using default: ${baseline_file}"
    }
    baseline_ch = channel.fromPath(baseline_file, checkIfExists: true)

    def ots28_baseline_file
    if (params.containsKey('ots28_baseline')) {
        ots28_baseline_file = params.ots28_baseline
    } else {
        ots28_baseline_file = "${projectDir}/data/baselines/full/RoSA_baseline_partial_infile.txt"
        log.info "No OTS28 baseline for Structure provided, using default: ${baseline_file}"
    }
    baseline_ch = channel.fromPath(baseline_file, checkIfExists: true)
    ots28_baseline_ch = channel.fromPath(ots28_baseline_file, checkIfExists: true)

    def rosa_allele_key_file
    if (params.containsKey('rosa_allele_key')) {
        rosa_allele_key_file = params.rosa_allele_key
    } else {
        rosa_allele_key_file = "${projectDir}/data/baselines/full/greb1_roha_alleles_reordered_wr.txt"
        log.info "No allele key provided, using default: ${rosa_allele_key_file}"
    }
    rosa_allele_key_ch = channel.fromPath(rosa_allele_key_file, checkIfExists: true)

    // Define input channels with automatic read type detection
    if (params.input_format == 'paired') {
        Channel
            .fromFilePairs(params.input, checkIfExists: true)
            .map { id, reads -> tuple(id, 'paired', reads) }
            .set { ch_input_reads }
    } else if (params.input_format == 'single') {
        Channel
            .fromPath(params.input, checkIfExists: true)
            .map { file -> tuple(file.simpleName, 'single', file) }
            .set { ch_input_reads }
    } else {
        // Auto-detect format based on file pattern matching
        Channel
            .fromFilePairs(params.input, checkIfExists: true, size: -1)
            .map { id, reads -> 
                def readType = reads.size() > 1 ? 'paired' : 'single'
                tuple(id, readType, reads)
            }
            .set { ch_input_reads }
    }
        
    ch_input_fastq = Channel
        .fromPath(params.input, checkIfExists: true)
        .collect()

        // Branch workflow based on read type
    ch_input_reads
        .branch {
            paired: it[1] == 'paired'
            single: it[1] == 'single'
        }
        .set { ch_reads_branched }

    // Process paired-end reads
    ch_reads_branched.paired
        .combine(adapters_ch)
        .set { ch_paired_adapters }


    
    FASTQC(ch_input_fastq) // FASTQC all input files
    TRIMMOMATIC(ch_paired_adapters, params.trim_params)
    FLASH2(TRIMMOMATIC.out.trimmed_paired, params.min_overlap, params.max_overlap)
    
    // Process single-end reads
    ch_reads_branched.single
        .combine(adapters_ch)
        .set { ch_single_adapters }
    
    TRIMMOMATIC_SINGLE(ch_single_adapters, params.trim_params)

    // Merge processed reads for downstream analysis
    ch_processed_reads = Channel.empty()
    ch_processed_reads = ch_processed_reads.mix(
        FLASH2.out.merged,
        TRIMMOMATIC_SINGLE.out.trimmed
    )
    

    bwa_input = ch_processed_reads.combine(reference_ch)
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

    // Group BAM and BAI files by reference
    bam_by_ref = SAMTOOLS.out.sorted_bam
        .map { sample_id, reference, bam -> 
            tuple(reference, bam) 
        }
        .groupTuple(by: 0)

    bai_by_ref = SAMTOOLS.out.sorted_bam_index
        .map { sample_id, reference, bai -> 
            tuple(reference, bai) 
        }
        .groupTuple(by: 0)

    // Join BAM and BAI groups by reference
    combined_bam_files = bam_by_ref
        .join(bai_by_ref)
    
    mpileup_input = combined_bam_files
        .map { reference, bams, bais -> 
            def ref_file = reference_files
                .collect { file(it) }
                .find { it.simpleName == reference } // Find the reference file by name. Index files have a simplename that does not match.
            if (!ref_file) {
                error "No exact match found for reference: ${reference}"
            }
            tuple(reference, bams, bais, ref_file)
        }

    BCFTOOLS_MPILEUP(mpileup_input)
    GREB_HAPSTR(BCFTOOLS_MPILEUP.out.filtered_vcf, rosa_allele_key_ch)

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
    //CHECK_FILE_UPDATE(HAP2GENO.out.new_index, locus_index_ch) | view

    // Run Rubias analyses
    RUN_RUBIAS(GREB_HAPSTR.out.ots28_report, HAP2GENO.out.numgeno, baseline_ch, params.panel.toLowerCase(), HAP2GENO.out.geno)
}
