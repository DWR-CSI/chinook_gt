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

// Validate required parameters
if (!params.outdir) {
    error "ERROR: --outdir parameter is required. Please specify an output directory."
}
if (!params.project) {
    error "ERROR: --project parameter is required. Please specify a project name."
}

// Set default values for thresholds if not specified in configs
params.ots28_missing_threshold = params.ots28_missing_threshold ?: 0.5
params.gsi_missing_threshold = params.gsi_missing_threshold ?: 0.6
params.pofz_threshold = params.pofz_threshold ?: 0.8
params.concat_all_reads = params.concat_all_reads ?: false
params.use_sequoia = params.use_sequoia ?: false
params.sequoia_mode = params.sequoia_mode ?: 'par'
params.sequoia_missing_threshold = params.sequoia_missing_threshold ?: 0.5
params.species_max_repro_age = params.species_max_repro_age ?: 7
params.species_min_repro_age = params.species_min_repro_age ?: 1
params.haplotype_depth = params.haplotype_depth ?: 4
params.total_depth = params.total_depth ?: 8
params.allele_balance = params.allele_balance ?: 0.35
params.loci_to_remove = params.loci_to_remove ?: ""
params.male_sexid_threshold = params.male_sexid_threshold ?: 0.02
params.female_sexid_threshold = params.female_sexid_threshold ?: 0.002
params.sexid_min_reads = params.sexid_min_reads ?: 10000
params.offspring_max_age = params.offspring_max_age ?: 6
params.offspring_maxBY = params.offspring_maxBY ?: params.offspring_birthyear
params.offspring_minBY = params.offspring_minBY ?: ((params.offspring_maxBY) ? (params.offspring_maxBY - params.offspring_max_age) : null)

// Validate numeric threshold ranges
if (params.ots28_missing_threshold < 0 || params.ots28_missing_threshold > 1) {
    error "ERROR: ots28_missing_threshold must be between 0 and 1 (got ${params.ots28_missing_threshold})"
}
if (params.gsi_missing_threshold < 0 || params.gsi_missing_threshold > 1) {
    error "ERROR: gsi_missing_threshold must be between 0 and 1 (got ${params.gsi_missing_threshold})"
}
if (params.pofz_threshold < 0 || params.pofz_threshold > 1) {
    error "ERROR: pofz_threshold must be between 0 and 1 (got ${params.pofz_threshold})"
}
if (params.allele_balance < 0 || params.allele_balance > 1) {
    error "ERROR: allele_balance must be between 0 and 1 (got ${params.allele_balance})"
}
if (params.male_sexid_threshold < 0 || params.male_sexid_threshold > 1) {
    error "ERROR: male_sexid_threshold must be between 0 and 1 (got ${params.male_sexid_threshold})"
}
if (params.female_sexid_threshold < 0 || params.female_sexid_threshold > 1) {
    error "ERROR: female_sexid_threshold must be between 0 and 1 (got ${params.female_sexid_threshold})"
}
if (params.haplotype_depth < 0) {
    error "ERROR: haplotype_depth must be non-negative (got ${params.haplotype_depth})"
}
if (params.total_depth < 0) {
    error "ERROR: total_depth must be non-negative (got ${params.total_depth})"
}
if (params.sexid_min_reads < 0) {
    error "ERROR: sexid_min_reads must be non-negative (got ${params.sexid_min_reads})"
}

if (params.use_sequoia) { // only validated if Sequoia is used
    if (params.sequoia_missing_threshold < 0 || params.sequoia_missing_threshold > 1) {
        error "ERROR: sequoia_missing_threshold must be between 0 and 1 (got ${params.sequoia_missing_threshold})"
    }
    if (params.species_max_repro_age < params.species_min_repro_age) {
        error "ERROR: species_max_repro_age (${params.species_max_repro_age}) must be >= species_min_repro_age (${params.species_min_repro_age})"
    }
    if (!params.offspring_maxBY) {
        error "ERROR: offspring_maxBY is required when use_sequoia is true"
    }
    if (params.offspring_maxBY && params.offspring_minBY && (params.offspring_maxBY < params.offspring_minBY)) {
        error "ERROR: offspring_maxBY (${params.offspring_maxBY}) must be >= offspring_minBY (${params.offspring_minBY})"
    }
}

// Import modules
include { FASTQC } from './modules/fastqc'
include { TRIMMOMATIC; TRIMMOMATIC_SINGLE } from './modules/trimmomatic'
include { FLASH2 } from './modules/flash2'
include { BWA_MEM } from './modules/bwa_mem'
include { SAMTOOLS } from './modules/samtools'
include { MULTIQC } from './modules/multiqc'
include { INDEX_REFERENCE } from './modules/index_reference'
include { ANALYZE_IDXSTATS } from './modules/idxstats_analysis'
include { GEN_MHP_SAMPLE_SHEET; PREP_MHP_RDS; GEN_HAPS} from './modules/microhaplot.nf'
include { RUN_RUBIAS } from './modules/rubias.nf'
include { STRUC_PARAMS; STRUCTURE } from './modules/structure.nf'
include { STRUCTURE_ROSA_REPORT } from './modules/rosa.nf'
include { BCFTOOLS_MPILEUP } from './modules/bcftools.nf'
include { GREB_HAPSTR } from './modules/RoSA_hap_str.nf'
include { CONCAT_READS } from './modules/concat_reads.nf'
include { RUN_SEQUOIA } from './modules/sequoia.nf'

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
    
    Microhaplotopia filtering parameters:
    Haplotype Depth Filter  : ${params.haplotype_depth}
    Haplotype Total Depth   : ${params.total_depth}
    Allele Balance Filter   : ${params.allele_balance}
    Rubias parameters:
    OTS28 Missing Threshold : ${params.ots28_missing_threshold}
    GSI Missing Threshold   : ${params.gsi_missing_threshold}
    PofZ Threshold          : ${params.pofz_threshold}
    Sex ID Parameters:
    Male Sex ID Threshold   : ${params.male_sexid_threshold}
    Female Sex ID Threshold : ${params.female_sexid_threshold}
    Sex ID Min Reads        : ${params.sexid_min_reads}
    Other parameters:
    Concatenate All Reads   : ${params.concat_all_reads}
    Sequoia Parameters:
    Use Sequoia             : ${params.use_sequoia}
    Sequoia Mode            : ${params.sequoia_mode}
    Sequoia Missing Threshold: ${params.sequoia_missing_threshold}
    Species Max Repro Age   : ${params.species_max_repro_age}
    Species Min Repro Age   : ${params.species_min_repro_age}
    Offspring Birth Year    : ${params.offspring_birthyear ?: 'Not specified'}
    Offspring Min Birth Year: ${params.offspring_minBY ?: 'Not specified'}
    Offspring Max Birth Year: ${params.offspring_maxBY ?: 'Not specified'}
    Loci Removal Regex    : ${params.loci_to_remove ?: 'None specified'}

    ==============================================
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
        log.info "No OTS28 baseline for Structure provided, using default: ${ots28_baseline_file}"
    }
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

    // Choose read processing approach based on parameter
    if (params.concat_all_reads) {
        // Collect all available reads per sample for concatenation
        ch_all_reads_per_sample = Channel.empty()
        
        // Add merged reads from FLASH2
        ch_all_reads_per_sample = ch_all_reads_per_sample.mix(
            FLASH2.out.merged
        )
        
        // Add unmerged reads from FLASH2 - transpose to separate paired files
        ch_all_reads_per_sample = ch_all_reads_per_sample.mix(
            FLASH2.out.unmerged.transpose()
        )
        
        // Add unpaired reads from TRIMMOMATIC - transpose to separate files
        ch_all_reads_per_sample = ch_all_reads_per_sample.mix(
            TRIMMOMATIC.out.trimmed_unpaired.transpose()
        )
        
        // Add single-end trimmed reads
        ch_all_reads_per_sample = ch_all_reads_per_sample.mix(
            TRIMMOMATIC_SINGLE.out.trimmed
        )
        
        // Group all files by sample_id and concatenate
        ch_grouped_reads = ch_all_reads_per_sample
            .groupTuple(by: 0)
        
        CONCAT_READS(ch_grouped_reads)
        ch_processed_reads = CONCAT_READS.out.concatenated
    } else {
        // Original approach - only merged and single-end reads
        ch_processed_reads = Channel.empty()
        ch_processed_reads = ch_processed_reads.mix(
            FLASH2.out.merged,
            TRIMMOMATIC_SINGLE.out.trimmed
        )
    }
    

    bwa_input = ch_processed_reads.combine(reference_ch)
    BWA_MEM(bwa_input)
    SAMTOOLS(BWA_MEM.out.aligned_sam)
    // Collect all idxstats files
    idxstats_files_sorted = SAMTOOLS.out.idxstats
        .branch { 
            main: it[1] ==~ /transition|full/
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
                error "No exact match found for reference: ${reference}. Use 'full' as the reference panel."
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

    // Run Rubias analyses
    RUN_RUBIAS(GREB_HAPSTR.out.ots28_report, baseline_ch, params.panel.toLowerCase(), GEN_HAPS.out.haps, ANALYZE_IDXSTATS.out.sexid)

    // Run PBT analysis
    if (params.use_sequoia) {
        par_geno_file = Channel.fromPath(params.parent_geno_input, checkIfExists: true)
        par_lh_file = Channel.fromPath(params.parent_lifehistory, checkIfExists: true)
        RUN_SEQUOIA(par_geno_file, par_lh_file, GEN_HAPS.out.haps, params.offspring_birthyear, params.offspring_minBY, params.offspring_maxBY, params.offspring_max_age)
    }
}
