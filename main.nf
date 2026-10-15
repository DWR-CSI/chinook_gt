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
include { DIMER_ANALYSIS } from './modules/dimer_analysis'
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
include { BCFTOOLS_MPILEUP; BCFTOOLS_MERGE } from './modules/bcftools.nf'
include { GREB_HAPSTR } from './modules/RoSA_hap_str.nf'
include { CONCAT_READS } from './modules/concat_reads.nf'
include { RUN_SEQUOIA } from './modules/sequoia.nf'
include { CKMR_PO } from './modules/CKMRsim.nf'
include { CKMRSIM_RUBIAS_SUMMARY } from './modules/summary.nf'
// Full genome mapping modules (for LFAR, WRAP, VGLL3SIX6 loci)
include { DOWNLOAD_AND_INDEX_GENOME } from './modules/fullgenome/download_genome'
include { MAKE_THINNED_GENOME; MAKE_THINNED_GENOME_MOUNT } from './modules/fullgenome/make_thinned_genome'
include { MAP_TO_FULL_GENOME; MAP_TO_FULL_GENOME_MOUNT } from './modules/fullgenome/map_fullgenome'
include { EXTRACT_READS_FROM_REGIONS } from './modules/fullgenome/extract_regions'
include { BAM_TO_FASTQ } from './modules/fullgenome/bam_to_fastq'
include { REMAP_TO_THINNED_GENOME } from './modules/fullgenome/remap_thinned'

// Helper functions for parameter type coercion (Nextflow 26+ string compatibility)
def toBool(val, defaultVal = false) {
    if (val == null) return defaultVal
    if (val instanceof Boolean) return val
    if (val instanceof String) {
        def s = val.trim().toLowerCase()
        if (s in ['true', '1', 'yes']) return true
        if (s in ['false', '0', 'no']) return false
    }
    return val ? true : false
}

def toNum(val, defaultVal = null) {
    if (val == null) return defaultVal
    if (val instanceof Number) return val
    if (val instanceof String) {
        try {
            return val.contains('.') ? Double.parseDouble(val.trim()) : Long.parseLong(val.trim())
        } catch (Exception _e) {
            return defaultVal
        }
    }
    return defaultVal
}

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
        def panel_lower = params.panel.toLowerCase()
        if (panel_lower == 'transition') {
                reference_files = [
                    "$projectDir/data/targets/transition/transition.fasta",
                    //"$projectDir/data/targets/LFAR/LFAR.fasta", // Commented out. Not used in DWR primers?
                    "$projectDir/data/targets/WRAP/WRAP.fasta"
                ]
                log.info "Using transition panel references"
        } else if (panel_lower == 'full') {
                reference_files = [
                    "$projectDir/data/targets/full/full.fna"//,
                    //"$projectDir/data/targets/full_VGLL3Six6LFARWRAP/VGLL3Six6LFARWRAP.fna"
                ]
                log.info "Using full panel reference"
        } else {
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


// Set default values for thresholds if not specified in configs


workflow {
// Check if user accidentally passed a parameter file with double-dash --params-file
if (params.containsKey('paramsFile') || params.containsKey('params-file') || params.containsKey('params_file')) {
    def pf = params.get('paramsFile') ?: params.get('params-file') ?: params.get('params_file')
    error """
    ================================================================================
    ERROR: Parameter file '${pf}' was specified with a double dash (--params-file).
    In Nextflow, core CLI options use a single dash: -params-file
    Please run:
      nextflow run main.nf -profile <profile> -params-file ${pf}
    ================================================================================
    """
}

// Validate required parameters
if (!params.input) {
    error "ERROR: --input parameter is required. Please specify an input file or pattern (e.g., via -params-file or --input 'path/*_R{1,2}.fastq.gz')."
}
if (!params.outdir) {
    error "ERROR: --outdir parameter is required. Please specify an output directory."
}
if (!params.project) {
    error "ERROR: --project parameter is required. Please specify a project name."
}

// Coerce booleans (handles Nextflow 26+ CLI string parsing)
def run_fastqc = toBool(params.run_fastqc)
def run_dimer_check = toBool(params.run_dimer_check, true)
def concat_all_reads = toBool(params.concat_all_reads)
def use_sequoia = toBool(params.use_sequoia)
def use_CKMR = toBool(params.use_CKMR)

// Coerce numbers for validation
def ots28_missing_threshold = toNum(params.ots28_missing_threshold)
def gsi_missing_threshold = toNum(params.gsi_missing_threshold)
def pofz_threshold = toNum(params.pofz_threshold)
def allele_balance = toNum(params.allele_balance)
def male_sexid_threshold = toNum(params.male_sexid_threshold)
def female_sexid_threshold = toNum(params.female_sexid_threshold)
def haplotype_depth = toNum(params.haplotype_depth)
def total_depth = toNum(params.total_depth)
def sexid_min_reads = toNum(params.sexid_min_reads)
def sequoia_missing_threshold = toNum(params.sequoia_missing_threshold)
def species_max_repro_age = toNum(params.species_max_repro_age)
def species_min_repro_age = toNum(params.species_min_repro_age)

// Validate numeric threshold ranges
if (ots28_missing_threshold != null && (ots28_missing_threshold < 0 || ots28_missing_threshold > 1)) {
    error "ERROR: ots28_missing_threshold must be between 0 and 1 (got ${params.ots28_missing_threshold})"
}
if (gsi_missing_threshold != null && (gsi_missing_threshold < 0 || gsi_missing_threshold > 1)) {
    error "ERROR: gsi_missing_threshold must be between 0 and 1 (got ${params.gsi_missing_threshold})"
}
if (pofz_threshold != null && (pofz_threshold < 0 || pofz_threshold > 1)) {
    error "ERROR: pofz_threshold must be between 0 and 1 (got ${params.pofz_threshold})"
}
if (allele_balance != null && (allele_balance < 0 || allele_balance > 1)) {
    error "ERROR: allele_balance must be between 0 and 1 (got ${params.allele_balance})"
}
if (male_sexid_threshold != null && (male_sexid_threshold < 0 || male_sexid_threshold > 1)) {
    error "ERROR: male_sexid_threshold must be between 0 and 1 (got ${params.male_sexid_threshold})"
}
if (female_sexid_threshold != null && (female_sexid_threshold < 0 || female_sexid_threshold > 1)) {
    error "ERROR: female_sexid_threshold must be between 0 and 1 (got ${params.female_sexid_threshold})"
}
if (haplotype_depth != null && haplotype_depth < 0) {
    error "ERROR: haplotype_depth must be non-negative (got ${params.haplotype_depth})"
}
if (total_depth != null && total_depth < 0) {
    error "ERROR: total_depth must be non-negative (got ${params.total_depth})"
}
if (sexid_min_reads != null && sexid_min_reads < 0) {
    error "ERROR: sexid_min_reads must be non-negative (got ${params.sexid_min_reads})"
}

if (use_sequoia) { // only validated if Sequoia is used
    if (sequoia_missing_threshold != null && (sequoia_missing_threshold < 0 || sequoia_missing_threshold > 1)) {
        error "ERROR: sequoia_missing_threshold must be between 0 and 1 (got ${params.sequoia_missing_threshold})"
    }
    if (species_max_repro_age != null && species_min_repro_age != null && species_max_repro_age < species_min_repro_age) {
        error "ERROR: species_max_repro_age (${params.species_max_repro_age}) must be >= species_min_repro_age (${params.species_min_repro_age})"
    }
    if (!params.offspring_maxBY) {
        error "ERROR: offspring_maxBY is required when use_sequoia is true"
    }
}
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
    Run FastQC              : ${params.run_fastqc}
    Run Dimer Check         : ${params.run_dimer_check}
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
        reference_ch = channel.empty()
    } else {
        channel
            .fromList(reference_files)
        .flatMap { ref ->
            // For each reference, check all required index files
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
            files.groupBy { f -> f.simpleName }
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
    _locus_index_ch = channel.fromPath(locus_index_file, checkIfExists: true)

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
    _ots28_baseline_ch = channel.fromPath(ots28_baseline_file, checkIfExists: true)

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
        channel
            .fromFilePairs(params.input, checkIfExists: true)
            .map { id, reads -> tuple(id, 'paired', reads) }
            .set { ch_input_reads }
    } else if (params.input_format == 'single') {
        channel
            .fromPath(params.input, checkIfExists: true)
            .map { file -> tuple(file.simpleName, 'single', file) }
            .set { ch_input_reads }
    } else {
        // Auto-detect format based on file pattern matching
        channel
            .fromFilePairs(params.input, checkIfExists: true, size: -1)
            .map { id, reads -> 
                def readType = reads.size() > 1 ? 'paired' : 'single'
                tuple(id, readType, reads)
            }
            .set { ch_input_reads }
    }
        
    if (run_fastqc) {
        ch_input_fastq = channel
            .fromPath(params.input, checkIfExists: true)
            .collect()
        FASTQC(ch_input_fastq) // FASTQC all input files
    }

    // Branch workflow based on read type
    ch_input_reads
        .branch { read ->
            paired: read[1] == 'paired'
            single: read[1] == 'single'
        }
        .set { ch_reads_branched }

    // Process paired-end reads
    ch_reads_branched.paired
        .combine(adapters_ch)
        .set { ch_paired_adapters }

    if (run_dimer_check) {
        DIMER_ANALYSIS(ch_reads_branched.paired, params.min_overlap, params.min_outie_overlap, params.max_overlap)

        merged_counts = DIMER_ANALYSIS.out.counts
            .collectFile(
                name: 'dimer_counts_mqc.tsv',
                keepHeader: true,
                storeDir: "${params.outdir}/${params.project}/dimers/summary"
            )
    } else {
        merged_counts = channel.empty()
    }

    TRIMMOMATIC(ch_paired_adapters, params.trim_params)
    FLASH2(TRIMMOMATIC.out.trimmed_paired, params.min_overlap, params.min_outie_overlap, params.max_overlap)
    
    // Process single-end reads
    ch_reads_branched.single
        .combine(adapters_ch)
        .set { ch_single_adapters }
    
    TRIMMOMATIC_SINGLE(ch_single_adapters, params.trim_params)

    // Choose read processing approach based on parameter
    if (concat_all_reads) {
        // Collect all available reads per sample for concatenation
        ch_all_reads_per_sample = channel.empty()
        
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
        ch_processed_reads = channel.empty()
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
    def ch_thinned_fasta = channel.empty()
    if (params.panel?.toLowerCase() == 'full') {
        log.info "Running parallel full genome mapping for VGLL3SIX6/LFAR/WRAP loci"
        log.info "Using combined region file: ${params.fullgenome_region_file}"

        // Initialize channels for merging downstream
        ch_thinned_fasta = channel.empty()
        def ch_map_bam = channel.empty()
        def ch_thinned_indices = channel.empty() // Will collect indices differently

        // Get the combined region file (contains all LFAR + WRAP + VGLL3SIX6 regions)
        region_file = channel.fromPath(params.fullgenome_region_file).collect()

        if (params.full_genome_mount_path) {
            log.info "Using mounted genome at: ${params.full_genome_mount_path}"
            
            // Create thinned genome using mounted files
            MAKE_THINNED_GENOME_MOUNT(
                params.fullgenome_ref_name,
                params.full_genome_mount_path,
                region_file
            )
            
            // Map merged reads using mounted files
            // Chunk inputs to reduce genome copy overhead
            ch_processed_reads_chunked = ch_processed_reads
                .toSortedList({ a, b -> a[0] <=> b[0] })
                .flatMap { sorted_list -> 
                    sorted_list.collate(params.fullgenome_chunk_size ?: 30).collect { chunk -> 
                        tuple(chunk.collect { c -> c[0] }, chunk.collect { c -> c[1] }) 
                    }
                }

            MAP_TO_FULL_GENOME_MOUNT(
                ch_processed_reads_chunked,
                params.full_genome_mount_path
            )

            ch_thinned_fasta = MAKE_THINNED_GENOME_MOUNT.out.fasta
            ch_map_bam = MAP_TO_FULL_GENOME_MOUNT.out.bam.transpose()

            // Collect indices from Mount version
            ch_thinned_indices = MAKE_THINNED_GENOME_MOUNT.out.indices

        } else {
            // Standard approach: Download and Index
            DOWNLOAD_AND_INDEX_GENOME()
            
            ch_genome_indices = DOWNLOAD_AND_INDEX_GENOME.out.indices

            // Create thinned genome
            MAKE_THINNED_GENOME(
                params.fullgenome_ref_name,
                DOWNLOAD_AND_INDEX_GENOME.out.genome,
                DOWNLOAD_AND_INDEX_GENOME.out.fai,
                region_file
            )

            // Map merged reads
            // Chunk inputs to reduce genome copy overhead
            ch_processed_reads_chunked = ch_processed_reads
                .toSortedList({ a, b -> a[0] <=> b[0] })
                .flatMap { sorted_list -> 
                    sorted_list.collate(params.fullgenome_chunk_size ?: 30).collect { chunk -> 
                        tuple(chunk.collect { c -> c[0] }, chunk.collect { c -> c[1] }) 
                    }
                }

            MAP_TO_FULL_GENOME(
                ch_processed_reads_chunked,
                DOWNLOAD_AND_INDEX_GENOME.out.genome,
                ch_genome_indices
            )

            ch_thinned_fasta = MAKE_THINNED_GENOME.out.fasta
            ch_map_bam = MAP_TO_FULL_GENOME.out.bam.transpose()
            
            // Collect indices from Standard version
            ch_thinned_indices = MAKE_THINNED_GENOME.out.indices
        }



        // Extract reads from target regions (works with bam from either path)
        EXTRACT_READS_FROM_REGIONS(
            ch_map_bam,
            params.fullgenome_ref_name,
            region_file
        )

        // Convert extracted BAM to FASTQ
        BAM_TO_FASTQ(EXTRACT_READS_FROM_REGIONS.out.bam)

        // Remap to thinned genome for correct coordinates
        REMAP_TO_THINNED_GENOME(
            BAM_TO_FASTQ.out.fastq,
            ch_thinned_fasta,
            ch_thinned_indices
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
        .branch { record ->
            main: record[1] ==~ "transition|full|${params.fullgenome_ref_name}"
            other: true
        }

    idxstats_main = idxstats_files_sorted.main
        .map { stat -> stat[2] }  // Extract just the file path from the tuple
        .collect()                // Collect all files into a single list

    // Run the analysis
    ANALYZE_IDXSTATS(idxstats_main)

    // Join sorted BAM and BAI by [sample_id, reference] to guarantee correct pairing
    SAMTOOLS.out.sorted_bam
        .join(SAMTOOLS.out.sorted_bam_index, by: [0, 1])
        .map { sample_id, reference, bam, bai -> 
            tuple(reference, sample_id, bam, bai) 
        }
        .groupTuple(by: 0)
        .flatMap { reference, sample_ids, bams, bais ->
            // Zip sample attributes and sort deterministically by sample_id to prevent cache invalidation on resume
            def paired = [sample_ids, bams, bais].transpose().sort { a, b -> a[0] <=> b[0] }
            def chunks = paired.collate(params.mpileup_chunk_size ?: 500)
            
            def result = []
            chunks.eachWithIndex { chunk, index ->
                def chunk_bams = chunk.collect { c -> c[1] }
                def chunk_bais = chunk.collect { c -> c[2] }
                result << tuple(reference, chunk_bams, chunk_bais, index)
            }
            return result
        }
        .set { mpileup_chunks }

    // Build mpileup_input - need to match reference files to BAM groups
    // For 'full' panel, we have two reference types: "full" and "Chinook_FullPanel_VGLL3SIX6LFARWRAP"
    if (params.panel?.toLowerCase() == 'full') {
        // Split by reference type
        mpileup_chunks
            .branch { chunk ->
                fullgenome: chunk[0] == params.fullgenome_ref_name
                targetfasta: true
            }
            .set { bam_branches }

        // Target FASTA BAMs - match to reference_files
        mpileup_targetfasta = bam_branches.targetfasta
            .map { reference, bams, bais, index ->
                def ref_file = reference_files
                    .collect { ref -> file(ref) }
                    .find { ref -> ref.simpleName == reference }
                if (!ref_file) {
                    error "No exact match found for reference: ${reference}"
                }
                tuple(reference, bams, bais, ref_file, index)
            }

        // Full genome BAMs - use thinned genome
        mpileup_fullgenome = bam_branches.fullgenome
            .combine(ch_thinned_fasta)
            .map { reference, bams, bais, index, ref_fasta ->
                tuple(reference, bams, bais, ref_fasta, index)
            }

        // Merge both
        mpileup_input = mpileup_targetfasta.mix(mpileup_fullgenome)
    } else {
        // Transition panel - only target FASTA references
        mpileup_input = mpileup_chunks
            .map { reference, bams, bais, index ->
                def ref_file = reference_files
                    .collect { ref -> file(ref) }
                    .find { ref -> ref.simpleName == reference }
                if (!ref_file) {
                    error "No exact match found for reference: ${reference}. Use 'full' as the reference panel."
                }
                tuple(reference, bams, bais, ref_file, index)
            }
    }

    BCFTOOLS_MPILEUP(mpileup_input)

    // Group MPILEUP outputs by reference and sort deterministically by filename before merging
    BCFTOOLS_MPILEUP.out.vcf
        .groupTuple(by: 0)
        .map { reference, vcfs, csis ->
            def sorted = [vcfs, csis].transpose().sort { a, b -> a[0].name <=> b[0].name }
            tuple(reference, sorted.collect { s -> s[0] }, sorted.collect { s -> s[1] })
        }
        .set { ch_vcfs_to_merge }

    BCFTOOLS_MPILEUP.out.filtered_vcf
        .groupTuple(by: 0)
        .map { reference, vcfs, csis ->
            def sorted = [vcfs, csis].transpose().sort { a, b -> a[0].name <=> b[0].name }
            tuple(reference, sorted.collect { s -> s[0] }, sorted.collect { s -> s[1] })
        }
        .set { ch_filtered_vcfs_to_merge }

    // Join them so they can be merged together
    ch_vcfs_to_merge
        .join(ch_filtered_vcfs_to_merge)
        .set { ch_merge_input }

    BCFTOOLS_MERGE(ch_merge_input)

    // Filter VCFs for Greb Hapstr (RoSA) - only need the main reference, not full genome
    BCFTOOLS_MERGE.out.filtered_vcf
        .filter { reference, _vcf -> 
             reference != params.fullgenome_ref_name 
        }
        .set { vcf_for_rosa }

    GREB_HAPSTR(vcf_for_rosa, rosa_allele_key_ch)

    // Group SAM files for microhaplotype analysis and chunk them
    ch_aligned_sam
        .groupTuple(by: 1)
        .flatMap { sample_ids, reference, sams ->
            // Pair each sample_id with its corresponding sam and sort deterministically
            def paired = [sample_ids, sams].transpose().sort { a, b -> a[0] <=> b[0] }
            // Chunk the paired list
            def chunks = paired.collate(params.mhp_chunk_size ?: 500)
            
            // Return a list of tuples for each chunk
            def result = []
            chunks.eachWithIndex { chunk, index ->
                def chunk_sample_ids = chunk.collect { c -> c[0] }
                def chunk_sams = chunk.collect { c -> c[1] }
                
                // Get the VCF file for this reference
                def vcfFile = file("${projectDir}/data/VCFs/${reference}.vcf")
                if (!vcfFile.exists()) {
                    error "VCF file not found for reference: ${reference}."
                }
                result << tuple(chunk_sample_ids, reference, vcfFile, chunk_sams, index)
            }
            return result
        }
        .set { packed_for_mhp_chunked }
    
    // Generate samplesheets for each chunk
    GEN_MHP_SAMPLE_SHEET(packed_for_mhp_chunked)

    // Run Microhaplot to generate RDS files for each chunk
    PREP_MHP_RDS(GEN_MHP_SAMPLE_SHEET.out.mhp_samplesheet)

    // Collect all RDS files across all chunks and references, sorting deterministically
    PREP_MHP_RDS.out.rds
        .map { _reference, rds_files -> rds_files }
        .collect()
        .map { all_rds_files -> 
             tuple("combined", all_rds_files.flatten().sort { a, b -> a.name <=> b.name }) 
        }
        .set { combined_rds }

    // Run GEN_HAPS on the combined set of RDS files
    GEN_HAPS(combined_rds)

    // Collect all QC files
    ch_multiqc_files = channel.empty()
    if (run_fastqc) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.fastqc_results)
    }
    ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log)
    ch_multiqc_files = ch_multiqc_files.mix(FLASH2.out.log)
    //ch_multiqc_files = ch_multiqc_files.mix(BWA_MEM.out.log)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS.out.idxstats.collect { s -> s[2] })
    if (run_dimer_check) {
        ch_multiqc_files = ch_multiqc_files.mix(merged_counts)
    }
    //ch_multiqc_files.view { println "Debug: MultiQC input file: $it" }
    MULTIQC(ch_multiqc_files.collect())
    
    RUN_RUBIAS(GREB_HAPSTR.out.ots28_report, baseline_ch, params.panel.toLowerCase(), GEN_HAPS.out.haps, ANALYZE_IDXSTATS.out.sexid)
    
    // Run CKMR analysis
    if (use_CKMR) {
        par_geno_wide = channel.fromPath(params.CKMR_parent_geno_input, checkIfExists: true)
        extra_genos_allele_freqs = channel.fromPath(params.CKMR_extra_genos_allele_freqs)
        CKMR_PO(GEN_HAPS.out.haps, par_geno_wide, extra_genos_allele_freqs)
        CKMRSIM_RUBIAS_SUMMARY(RUN_RUBIAS.out.summary, CKMR_PO.out.PO_results)
    }
    // Run PBT analysis
    if (use_sequoia) {
        par_geno_file = channel.fromPath(params.parent_geno_input, checkIfExists: true)
        par_lh_file = channel.fromPath(params.parent_lifehistory, checkIfExists: true)
        RUN_SEQUOIA(par_geno_file, par_lh_file, GEN_HAPS.out.haps, params.offspring_birthyear, params.offspring_minBY, params.offspring_maxBY, params.offspring_max_age)
    }
}
