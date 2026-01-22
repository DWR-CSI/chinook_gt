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
params.fullgenome_ref_name = params.fullgenome_ref_name ?: 'Chinook_FullPanel_VGLL3Six6LFARWRAP'
params.ots28_missing_threshold = params.ots28_missing_threshold ?: 0.5
params.full_genome_mount_path = params.full_genome_mount_path ?: null
params.gsi_missing_threshold = params.gsi_missing_threshold ?: 0.6
params.pofz_threshold = params.pofz_threshold ?: 0.8
params.concat_all_reads = params.concat_all_reads ?: false
params.use_sequoia = params.use_sequoia ?: false
params.sequoia_mode = params.sequoia_mode ?: 'par'
params.sequoia_missing_threshold = params.sequoia_missing_threshold ?: 0.5
params.species_max_repro_age = params.species_max_repro_age ?: 6
params.species_min_repro_age = params.species_min_repro_age ?: 1
params.haplotype_depth = params.haplotype_depth ?: 4
params.total_depth = params.total_depth ?: 8
params.allele_balance = params.allele_balance ?: 0.35
params.loci_to_remove = params.loci_to_remove ?: ""
params.male_sexid_threshold = params.male_sexid_threshold ?: 0.02
params.female_sexid_threshold = params.female_sexid_threshold ?: 0.002
params.sexid_min_reads = params.sexid_min_reads ?: 10000
params.offspring_max_age = params.offspring_max_age ?: 6
params.offspring_birthyear = params.offspring_birthyear ?: 'unknown'
params.offspring_maxBY = params.offspring_maxBY ?: params.offspring_birthyear
params.offspring_minBY = params.offspring_minBY ?: ((params.offspring_maxBY) ? (params.offspring_maxBY - params.offspring_max_age) : null)
params.use_CKMR = params.use_CKMR ?: false
params.CKMR_logl_threshold = params.CKMR_logl_threshold ?: 6.9
params.CKMR_min_loci = params.CKMR_min_loci ?: 90
params.CKMR_parent_geno_input = params.CKMR_parent_geno_input ?: "$projectDir/examples/PBT/FRH2024_reference_genotypes.csv"
params.CKMR_extra_genos_allele_freqs = params.CKMR_extra_genos_allele_freqs ?: "$projectDir/examples/PBT/JPE2022-2024_geno_wide.csv"
params.fullgenome_region_file = params.fullgenome_region_file ?: "$projectDir/data/regions/Chinook_FullPanel_VGLL3Six6LFARWRAP-Otsh_v1.0.txt"

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
include { CKMR_PO } from './modules/CKMRsim.nf'
include { CKMRSIM_RUBIAS_SUMMARY } from './modules/summary.nf'
// Full genome mapping modules (for LFAR, WRAP, VGLL3SIX6 loci)
include { DOWNLOAD_AND_INDEX_GENOME } from './modules/fullgenome/download_genome'
include { MAKE_THINNED_GENOME } from './modules/fullgenome/make_thinned_genome'
include { MAP_TO_FULL_GENOME } from './modules/fullgenome/map_fullgenome'
include { EXTRACT_READS_FROM_REGIONS } from './modules/fullgenome/extract_regions'
include { BAM_TO_FASTQ } from './modules/fullgenome/bam_to_fastq'
include { REMAP_TO_THINNED_GENOME } from './modules/fullgenome/remap_thinned'

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
    Full Genome Ref : ${params.fullgenome_ref_name ?: 'N/A'}
    Region File     : ${params.fullgenome_region_file ?: 'N/A'}
    Adapter File    : ${params.adapter_file ?: 'Default'}
    
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
    if (reference_files.isEmpty()) {
        reference_ch = Channel.empty()
    } else {
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
                // add an INDEX_REFERENCE process here in the future?
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
    }
    
    // Define input channels
    def adapter_file_final = params.adapter_file ?: "${projectDir}/data/adapters/GTseq-PE.fa"
    adapters_ch = channel.fromPath(adapter_file_final)

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

    // Standard target FASTA mapping workflow (always runs for transition and full panels)
    bwa_input = ch_processed_reads.combine(reference_ch)
    BWA_MEM(bwa_input)

    // For 'full' panel, also run full genome mapping for VGLL3SIX6/LFAR/WRAP loci
    // These loci require full genome mapping to correctly handle off-target amplicons
    if (params.panel?.toLowerCase() == 'full') {
        log.info "Running parallel full genome mapping for VGLL3SIX6/LFAR/WRAP loci"
        log.info "Using combined region file: ${params.fullgenome_region_file}"

        // Determine if using mount path or staged files
        def mount_path = params.full_genome_mount_path
        def use_mount = mount_path != null
        
        // Initialize channel variables
        def ch_genome = Channel.empty()
        def ch_genome_fai = Channel.empty()
        def ch_genome_indices = Channel.empty()

        if (use_mount) {
            log.info "Using mounted genome at: ${mount_path}"
            // Skip DOWNLOAD_AND_INDEX_GENOME
        } else {
            // Download and index the full genome (cached via storeDir)
            DOWNLOAD_AND_INDEX_GENOME()
            ch_genome = DOWNLOAD_AND_INDEX_GENOME.out.genome
            ch_genome_fai = DOWNLOAD_AND_INDEX_GENOME.out.fai
            ch_genome_indices = DOWNLOAD_AND_INDEX_GENOME.out.amb
                .mix(DOWNLOAD_AND_INDEX_GENOME.out.ann)
                .mix(DOWNLOAD_AND_INDEX_GENOME.out.bwt)
                .mix(DOWNLOAD_AND_INDEX_GENOME.out.pac)
                .mix(DOWNLOAD_AND_INDEX_GENOME.out.sa)
                .collect()
        }

        // Get the combined region file (contains all LFAR + WRAP + VGLL3SIX6 regions)
        region_file = Channel.fromPath(params.fullgenome_region_file).collect()

        // Create thinned genome using combined regions (cached via storeDir)
        // If use_mount is true, ch_genome and ch_genome_fai are empty, but mount_path is used
        MAKE_THINNED_GENOME(
            params.fullgenome_ref_name,
            ch_genome,
            ch_genome_fai,
            region_file
        )

        // Map merged reads to full genome
        MAP_TO_FULL_GENOME(
            ch_processed_reads,
            DOWNLOAD_AND_INDEX_GENOME.out.genome.collect(),
            genome_indices
        )

        // Extract reads from target regions
        EXTRACT_READS_FROM_REGIONS(
            MAP_TO_FULL_GENOME.out.bam,
            params.fullgenome_ref_name,
            region_file
        )

        // Convert extracted BAM to FASTQ
        BAM_TO_FASTQ(EXTRACT_READS_FROM_REGIONS.out.bam)

        // Collect thinned genome index files for remapping
        thinned_indices = MAKE_THINNED_GENOME.out.amb
            .mix(MAKE_THINNED_GENOME.out.ann)
            .mix(MAKE_THINNED_GENOME.out.bwt)
            .mix(MAKE_THINNED_GENOME.out.pac)
            .mix(MAKE_THINNED_GENOME.out.sa)
            .collect()

        // Remap to thinned genome for correct coordinates
        REMAP_TO_THINNED_GENOME(
            BAM_TO_FASTQ.out.fastq,
            MAKE_THINNED_GENOME.out.fasta.collect(),
            thinned_indices
        )

        // Merge SAM outputs from both workflows
        ch_aligned_sam = BWA_MEM.out.aligned_sam.mix(REMAP_TO_THINNED_GENOME.out.aligned_sam)
    } else {
        // Transition panel - only target FASTA mapping
        ch_aligned_sam = BWA_MEM.out.aligned_sam
    }

    SAMTOOLS(ch_aligned_sam)
    // Collect all idxstats files
    idxstats_files_sorted = SAMTOOLS.out.idxstats
        .branch {
            main: it[1] ==~ /transition|full|${params.fullgenome_ref_name}/
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

    // Build mpileup_input - need to match reference files to BAM groups
    // For 'full' panel, we have two reference types: "full" and "Chinook_FullPanel_VGLL3Six6LFARWRAP"
    if (params.panel?.toLowerCase() == 'full') {
        // Split by reference type
        combined_bam_files
            .branch {
                fullgenome: it[0] == params.fullgenome_ref_name
                targetfasta: true
            }
            .set { bam_branches }

        // Target FASTA BAMs - match to reference_files
        mpileup_targetfasta = bam_branches.targetfasta
            .map { reference, bams, bais ->
                def ref_file = reference_files
                    .collect { file(it) }
                    .find { it.simpleName == reference }
                if (!ref_file) {
                    error "No exact match found for reference: ${reference}"
                }
                tuple(reference, bams, bais, ref_file)
            }

        // Full genome BAMs - use thinned genome
        mpileup_fullgenome = bam_branches.fullgenome
            .combine(MAKE_THINNED_GENOME.out.fasta)
            .map { reference, bams, bais, ref_fasta ->
                tuple(reference, bams, bais, ref_fasta)
            }

        // Merge both
        mpileup_input = mpileup_targetfasta.mix(mpileup_fullgenome)
    } else {
        // Transition panel - only target FASTA references
        mpileup_input = combined_bam_files
            .map { reference, bams, bais ->
                def ref_file = reference_files
                    .collect { file(it) }
                    .find { it.simpleName == reference }
                if (!ref_file) {
                    error "No exact match found for reference: ${reference}. Use 'full' as the reference panel."
                }
                tuple(reference, bams, bais, ref_file)
            }
    }

    BCFTOOLS_MPILEUP(mpileup_input)

    // Filter VCFs for Greb Hapstr (RoSA) - only need the main reference, not full genome
    BCFTOOLS_MPILEUP.out.filtered_vcf
        .filter { reference, vcf -> 
             reference != params.fullgenome_ref_name 
        }
        .set { vcf_for_rosa }

    GREB_HAPSTR(vcf_for_rosa, rosa_allele_key_ch)

    // Group SAM files for microhaplotype analysis
    ch_aligned_sam
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
        .set { packed_for_mhp }
    
    // Generate samplesheets
    GEN_MHP_SAMPLE_SHEET(packed_for_mhp)

    // Run Microhaplot to generate RDS files
    PREP_MHP_RDS(GEN_MHP_SAMPLE_SHEET.out.mhp_samplesheet)

    // If we have multiple RDS files (e.g. from target and full genome panels), we need to merge them.
    // Collect all RDS files across all references
    PREP_MHP_RDS.out.rds
        .map { reference, rds_files -> rds_files }
        .collect()
        .map { all_rds_files -> 
             // Use a generic reference name for the combined set, or just use the project name
             tuple("combined", all_rds_files) 
        }
        .set { combined_rds }

    // Run GEN_HAPS on the combined set of RDS files
    GEN_HAPS(combined_rds)

    // Collect all QC files
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.fastqc_results)
    ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log)
    ch_multiqc_files = ch_multiqc_files.mix(FLASH2.out.log)
    //ch_multiqc_files = ch_multiqc_files.mix(BWA_MEM.out.log)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS.out.idxstats.collect{it[2]})
    //ch_multiqc_files.view { println "Debug: MultiQC input file: $it" }
    MULTIQC(ch_multiqc_files.collect())
    
    RUN_RUBIAS(GREB_HAPSTR.out.ots28_report, baseline_ch, params.panel.toLowerCase(), GEN_HAPS.out.haps, ANALYZE_IDXSTATS.out.sexid)
    
    // Run CKMR analysis
    if (params.use_CKMR) {
        par_geno_wide = Channel.fromPath(params.CKMR_parent_geno_input, checkIfExists: true)
        extra_genos_allele_freqs = Channel.fromPath(params.CKMR_extra_genos_allele_freqs)
        CKMR_PO(GEN_HAPS.out.haps, par_geno_wide, extra_genos_allele_freqs)
        CKMRSIM_RUBIAS_SUMMARY(RUN_RUBIAS.out.summary, CKMR_PO.out.PO_results)
    }
    // Run PBT analysis
    if (params.use_sequoia) {
        par_geno_file = Channel.fromPath(params.parent_geno_input, checkIfExists: true)
        par_lh_file = Channel.fromPath(params.parent_lifehistory, checkIfExists: true)
        RUN_SEQUOIA(par_geno_file, par_lh_file, GEN_HAPS.out.haps, params.offspring_birthyear, params.offspring_minBY, params.offspring_maxBY, params.offspring_max_age)
    }
}
