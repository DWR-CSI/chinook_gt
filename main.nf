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
//include { FLASH2 } from './modules/flash2'
//include { BWA_MEM } from './modules/bwa_mem'
//include { SAMTOOLS } from './modules/samtools'
include { MULTIQC } from './modules/multiqc'


// Define input channels
Channel // paired
    .fromFilePairs(params.input, checkIfExists: true)
    .set { ch_input_fastq_pairs }
ch_input_fastq = Channel // all fastq files for FastQC
    .fromPath(params.input, checkIfExists: true)
    .collect()


    
workflow {
  //      ch_reads = Channel
    //    .fromFilePairs(params.input, checkIfExists: true)
      //  .view { sample_id, files -> 
        //    "Sample: $sample_id, Files: $files, Paths: ${files.collect { it.toString() }}"
        //}
    FASTQC(ch_input_fastq)
    TRIMMOMATIC(ch_input_fastq_pairs, params.adapter_file, params.trim_params)
    //FLASH2(TRIMMOMATIC.out.trimmed_paired, params.min_overlap, params.max_overlap)
    //BWA_MEM(FLASH2.out.merged, params.reference)
    //SAMTOOLS(BWA_MEM.out.aligned_sam)
    
    // Collect all QC files
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.fastqc_results)
    //ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log)
    //ch_multiqc_files = ch_multiqc_files.mix(FLASH2.out.log)
    //ch_multiqc_files = ch_multiqc_files.mix(BWA_MEM.out.log)
    //ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS.out.idxstats)
    MULTIQC(ch_multiqc_files.collect())
}

