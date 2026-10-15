# chinook_gt

## Table of Contents
- [Overview](#overview)
- [Setup](#setup)
- [Usage](#usage)
- [Execution Profiles](#execution-profiles)
- [Parameters](#parameters)
  - [Input & Output Options](#input--output-options)
  - [Trimmomatic Options](#trimmomatic-options)
  - [FLASH2 Options](#flash2-options)
  - [Processing Options](#processing-options)
  - [Reference & Baseline Options](#reference--baseline-options)
  - [Full Genome Mapping Options](#full-genome-mapping-options)
  - [Microhaplotopia Options](#microhaplotopia-options)
  - [Threshold Options](#threshold-options)
  - [Sex ID Options](#sex-id-options)
  - [CKMRsim Options](#ckmrsim-options)
  - [Sequoia Options](#sequoia-options)
  - [Loci Removal](#loci-removal)
- [Output](#output)
- [Troubleshooting](#troubleshooting)
- [Contributors](#contributors)
- [License](#license)

## Overview

This Nextflow pipeline takes .FASTQ files as input and produces Chinook run identifications as output along with relevant intermediate files.

## Setup
Currently, this pipeline is optimized for use with Microsoft Azure Batch (and other clouds) but can/will be adapted for use with HPC systems by specifying Singularity containers instead of Docker containers.

Prior to running the primary `main.nf` workflow, simply ensure that the parameters and paths specified in the configuration files located at `nextflow.config`, `conf/azure.config` (or your preferred profile), and `~/.nextflow/config` are correct for your setup. You can then run `main.nf` with Nextflow while specifying a parameters file (`.yml` or `.json`) corresponding to your run. Fully documented example parameter files are provided in [examples/settings/test_settings.yml](examples/settings/test_settings.yml) and [examples/settings/test_settings_azure.yml](examples/settings/test_settings_azure.yml), with quick-start configs at [test_settings.yml](test_settings.yml), [test_settings_azure.yml](test_settings_azure.yml), and [test_settings.json](test_settings.json).

## Usage

Run the pipeline with Azure Batch:
```bash
nextflow main.nf -profile azure -params-file settings.yml
```

Run locally with Docker:
```bash
nextflow main.nf -profile docker -params-file settings.yml
```

Run with Singularity (for HPC systems):
```bash
nextflow main.nf -profile singularity -params-file settings.yml
```

## Execution Profiles

The pipeline supports multiple execution profiles configured in the `conf/` directory:

- **`azure`** - Azure Batch execution with Docker containers (default cloud profile)
  - Uses Azure Batch for job orchestration
  - Auto-scaling compute pools
  - Azure Blob Storage for work and results

- **`azure_large`** - Azure Batch profile optimized for large sample numbers (>1,000 samples)
  - Increased memory allocations for `GEN_HAPS`, `PREP_MHP_RDS`, `RUN_RUBIAS`, `MULTIQC`, and `ANALYZE_IDXSTATS` to prevent OOM errors

- **`azure_local`** - Azure Batch with local adapter files
  - Similar to `azure` but uses local reference files instead of cloud storage

- **`docker`** - Local execution with Docker containers
  - For running on local machines or servers with Docker installed

- **`ucdavis_hive`** - Pre-configured profile for UC Davis Hive HPC
  - Uses SLURM executor with Singularity
  - Optimized for the UC Davis Hive cluster environment

### Creating Custom Profiles

To create a custom profile for your environment:

1. Create a new `.config` file in the `conf/` directory (e.g., `conf/my_hpc.config`)
2. Define your executor, container engine, and resource settings:
   ```groovy
   singularity.enabled = true
   workDir = '/path/to/work/dir'

   process {
       executor = 'slurm'
       queue = 'your-queue-name'
       // Add process-specific resource allocations
   }
   ```
3. Reference your profile when running: `nextflow main.nf -c conf/my_hpc.config -params-file settings.yml`

## Parameters

### Input & Output Options

- **`project`** (string, default: `"YourProject1"`)
  Name of the project or analysis run. Must contain only alphanumeric characters, underscores, or hyphens (pattern: `^[A-Za-z0-9_-]+$`). This name is used as the output subdirectory under `outdir`.

- **`outdir`** (path / string, required)
  Path to the root output directory where all pipeline outputs and intermediate files will be stored. Can be a local filesystem directory (e.g., `'test/output'`) or an Azure Blob Storage container path (e.g., `'az://nextflow/results/chinookGT_full/'`).

- **`panel`** (string, default: `"full"`, options: `"full"` or `"transition"`)
  GT-seq panel version and reference strategy:
  - `"full"`: 204-locus panel. Employs target FASTA mapping alongside hybrid full-genome alignment for off-target loci (`VGLL3SIX6`, `LFAR`, `WRAP`).
  - `"transition"`: Earlier transition panel with fewer loci.

- **`n_loci`** (integer, default: `204`)
  Total number of expected loci in the panel (full panel has 204 loci). Used in `ANALYZE_IDXSTATS` for read distribution metrics and quality control calculations.

- **`input`** (path / string, default: `"az://seqera1/raw_data/ADPH_mushed/*_R{1,2}_concatenated.fastq.gz"`)
  Path or glob pattern pointing to raw paired-end (or single-end) FASTQ files. Supports gzip compression (`.fastq.gz`). For paired-end data, use `{1,2}` to specify forward and reverse read pairs. Can be a local path or cloud storage URI (`az://...`).

- **`input_format`** (string, default: `"auto"`, options: `"auto"`, `"paired"`, `"single"`)
  Input read format detection mode:
  - `"auto"`: Automatically detects paired-end vs. single-end reads based on pattern matching (recommended).
  - `"paired"`: Strictly treats matched files as paired-end reads.
  - `"single"`: Treats input files as single-end reads.

### Trimmomatic Options

- **`adapter_file`** (path, default: `"data/adapters/GTseq-PE.fa"`)
  Path to FASTA file containing sequencing adapter and primer sequences for Trimmomatic clipping. Staged automatically in containerized runs. Set to `null` to disable adapter trimming.

- **`trim_params`** (string, default: `"LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36"`)
  Trimmomatic quality trimming and filtering parameters:
  - `LEADING:15`: Remove low-quality bases from start (below Phred Q15).
  - `TRAILING:15`: Remove low-quality bases from end (below Phred Q15).
  - `SLIDINGWINDOW:4:15`: 4-base wide sliding window; cut when average quality drops below Q15.
  - `MINLEN:36`: Drop reads that fall below 36 bp after trimming.

### FLASH2 Options

- **`min_overlap`** (integer, default: `3`, minimum: 1)
  Minimum overlap length (bp) required for FLASH2 to merge forward and reverse paired-end reads. The default of 3 accommodates variable GT-seq amplicon lengths.

- **`max_overlap`** (integer, default: `400`, minimum: 1)
  Maximum overlap length (bp) for FLASH2 to merge paired-end reads.

- **`min_outie_overlap`** (integer, default: `20`, minimum: 1)
  Minimum overlap length (bp) for merging reads in "outie" orientation.

### Processing Options

- **`run_fastqc`** (boolean, default: `false`)
  Run FastQC quality control analysis on raw FASTQ files. Set to `true` to generate per-sample FastQC reports with MultiQC integration.

- **`run_dimer_check`** (boolean, default: `true`)
  Run primer-dimer analysis via the `DIMER_ANALYSIS` module. Set to `false` to disable. Recommended `true` for GT-seq amplicon pools.

- **`concat_all_reads`** (boolean, default: `false`)
  Concatenate all available reads per sample (merged, unmerged, and unpaired) for BWA mapping. If `false` (default), uses only merged and single-end reads.

- **`loci_to_remove`** (string, default: `""`)
  Perl-compatible regular expression to match loci column names to remove from RUBIAS and Sequoia analyses. See [Loci Removal](#loci-removal) section for details.

### Reference & Baseline Options

- **`locus_index`** (path, default: `"data/indices/transition_loci_indices_20241104.csv"`)
  Path to CSV file mapping panel loci names to convert alleles to numbers based on historical CDFW usage. Used for Structure analysis only (not recommended).

- **`baseline`** (path, default: `"data/baselines/full/SWFSC-chinook-reference-baseline-CV.csv"`)
  Reference baseline CSV file for Genetic Stock Identification (GSI) with RUBIAS. Contains reference allele frequencies across Central Valley Chinook salmon reporting units. Can be a local filesystem path or Azure Blob Storage path (`az://...`).

- **`ots28_baseline`** (path, default: `"data/baselines/full/RoSA_baseline_partial_infile.txt"`)
  Reference baseline file for Structure/RoSA (OTS28) run-timing and ecotype assignment (e.g., spring vs. fall run), if Structure is used (not recommended).

- **`rubias_show_missing_data`** (string / boolean, default: `"TRUE"`)
  Whether RUBIAS output should include missing data (not drop loci (`"TRUE"` or `"FALSE"`). Recommended.

### Full Genome Mapping Options

- **`genome_cache`** (path, default: `"${projectDir}/data/genome_cache"`)
  Path to the genome cache directory. Used to store the reference genome and BWA indices after they are downloaded and indexed for the first time. If this is blank, a default path will be used. If empty, genome will be downloaded and indexed.

- **`full_genome_mount_path`** (path, default: `null`)
  (Work in progress) (Not recommended) Used when running the full panel in large cohorts with Azure File Share mounts. This should **not** be set when not using Azure File Share.
  Specifies the absolute path to the mounted Azure File Share directory containing the reference genome and BWA indices.
  Example: `/mnt/batch/tasks/fsmounts/seqres/reference/genomes`


### Microhaplotopia Options

- **`haplotype_depth`** (integer, default: `4`)
  Remove haplotypes with depth less than this value to filter potential genotyping errors or alleles with questionable depths.

- **`total_depth`** (integer, default: `8`)
  Minimum total number of reads required for a locus genotype to be retained.

- **`allele_balance`** (number, default: `0.35`, range: 0-1)
  Minimum allele balance ratio required for heterozygote genotypes to be called.

### Threshold Options

- **`ots28_missing_threshold`** (number, default: `0.7`, range: 0-1)
  Missing data threshold for OTS28 (RoSA) run-timing analysis. If less than this proportion of OTS28 data is missing, consider OTS28 data Intermediate instead of uncertain (multiplied by 100 in analysis).

- **`gsi_missing_threshold`** (number, default: `0.6`, range: 0-1)
  Genetic Stock Identification (GSI) missing data threshold. If an individual has missing genotype data at more than this proportion of GSI loci, its stock assignment is considered invalid.

- **`pofz_threshold`** (number, default: `0.8`, range: 0-1)
  RUBIAS maximum posterior probability (Probability of Z-score) threshold. If the maximum PofZ is less than this value, the stock assignment is considered ambiguous.

### Sex ID Options

- **`male_sexid_threshold`** (number, default: `0.01`, range: 0-1)
  Minimum proportion of reads mapped to the Y-chromosome locus (`sdy_I183`) relative to total panel reads to classify a sample as male.

- **`female_sexid_threshold`** (number, default: `0.002`, range: 0-1)
  Maximum proportion of reads mapped to the Y-chromosome locus (`sdy_I183`) relative to total panel reads to classify a sample as female. Samples between female and male thresholds are classified as unknown.

- **`sexid_min_reads`** (number, default: `10000`, minimum: 0)
  Minimum total reads mapped across the panel required to assign a confident sex identification. Samples with fewer reads are classified as unknown sex.

### CKMRsim Options

- **`use_CKMR`** (boolean, default: `false`)
  Enable CKMRsim for close-kin mark-recapture parent-offspring analysis. When enabled, the pipeline runs CKMRsim after genotype generation and combines results with RUBIAS stock assignments.

- **`CKMR_parent_geno_input`** (path, default: `"$projectDir/examples/PBT/FRH2024_reference_genotypes.csv"`)
  Path to parent reference genotype file in wide-format CSV. File format requirements:
  Header row with `SAMPLE_ID` as first column, followed by locus names with `.1` and `.2` suffixes (e.g., `Locus1.1`, `Locus1.2`). Missing data should be coded as `NA` or blank.

- **`CKMR_logl_threshold`** (number, default: `6.9`, minimum: 0)
  Log-likelihood ratio threshold for identifying parent-offspring relationships. Results with log-likelihood ratios below this threshold are filtered out. Higher values are more conservative.

- **`CKMR_min_loci`** (integer, default: `90`, minimum: 1)
  Minimum number of genotyped loci required for a valid parent-offspring assignment.

- **`CKMR_extra_genos_allele_freqs`** (path, default: `"$projectDir/examples/PBT/JPE2022-2024_geno_wide.csv"`)
  Path to extra genotypes file (wide-format CSV) containing additional individuals to include in allele frequency calculations. Including a larger and more representative sample improves the accuracy of kinship likelihood calculations. Follows the same format as the parent genotype file.

### Sequoia Options

- **`use_sequoia`** (boolean, default: `false`)
  Enable Sequoia for parentage-based tagging (PBT) and pedigree analysis.

- **`sequoia_mode`** (string, default: `"par"`, options: `"par"` or `"ped"`)
  Sequoia analysis mode: `"par"` for parentage analysis only (faster, recommended for large cohorts) or `"ped"` for full pedigree analysis including siblings (slower, computationally intensive).

- **`sequoia_missing_threshold`** (number, default: `0.5`, range: 0-1)
  Maximum allowable proportion of missing loci for Sequoia. If more than this proportion of data is missing, consider the result invalid. Set to `0.2` for conservative cutoff.

- **`parent_geno_input`** (path, default: `"$projectDir/examples/PBT/FRH2024_reference_genotypes.csv"`)
  Path to reference parent genotype file formatted for Sequoia.

- **`parent_lifehistory`** (path, default: `"$projectDir/examples/PBT/FRH2024_reference_lifehistory.txt"`)
  Path to reference parent metadata/life-history file formatted for Sequoia (`ID`, `Sex`, `BirthYear`).

- **`offspring_birthyear`** (string / integer, default: `"unknown"`)
  Birth year of the offspring cohort. Set to an integer year or `"unknown"`.

- **`offspring_minBY`** (string / integer, default: `"unknown"`)
  Minimum birth year for offspring cohort.

- **`offspring_maxBY`** (string / integer, default: `"2025"`)
  Maximum birth year for offspring cohort. Set to sampling year or current year if unknown.

- **`offspring_max_age`** (integer, default: `8`)
  Maximum age (in years) of offspring in the dataset.

- **`species_min_repro_age`** (integer, default: `1`)
  Minimum reproductive age of the species (Chinook salmon default: 1).

- **`species_max_repro_age`** (integer, default: `8`)
  Maximum reproductive age of the species (Chinook salmon default: 8).

### Loci Removal

The pipeline supports removing specific loci from analysis using the `loci_to_remove` parameter. This uses Perl-compatible regular expressions to match loci column names.

**How it works:**
- The regex pattern is matched against these column names in both RUBIAS and Sequoia analyses
- Matched loci are removed before Rubias and Sequoia analyses begins, but after the point of haplotype/genotype generation.

**Configuration:**
In your parameters YAML file:
```yaml
loci_to_remove: "(NC_037130\\.1:1062935[.-]1063235)|(NC_037130\\.1:864908[.-]865208)|(NC_037104\\.1:56552952[.-]56553042)"
```

Note that, above, the pattern `[.-]` is used to match either a period or a dash since the rubias.R script converts dashes to periods, but the Sequoia script and others may not.

**Notes:**
- Leave empty (`""`) to remove no loci
- Use Perl-compatible regex syntax
- Cannot contain single quotes
- Backslashes must be escaped (use `\\.` to match a literal period)
- Invalid patterns will be caught and no loci will be removed

## Output

The pipeline produces various intermediate and final files in the output directory specified in your configuration. The output is organized into subdirectories by analysis type.

### Output Directory Structure

```
results (user-specified name)/
├── {project_name}/
│   ├── fastqc/              # Quality control reports for raw reads
│   ├── trimmomatic/         # Trimmed read files
│   ├── flash2/              # Merged paired-end reads
│   ├── bwa/                 # Alignment files (BAM)
│   ├── microhaplotopia/     # Microhaplotype calling outputs
│   ├── rubias/              # Genetic stock assignment results
│   ├── ckmrsim/             # Parent-offspring analysis with CKMRsim (if enabled)
│   ├── summary/             # Summary report files (if CKMRsim enabled)
│   ├── sequoia/             # Parentage analysis (if enabled)
│   └── multiqc/             # Aggregated QC report
```

### Key Output Files

**Primary Results (rubias/ directory):**
- **`*_summary.tsv`** - **Primary output file** containing final genetic run IDs with RoSA-informed assignments. Use this file for DWR Central Valley Chinook salmon run identification.
- `*_filtered_haplotype.csv` - Processed genotype data used for RUBIAS analysis
- `*_observed_unfiltered_haplotype.csv` - Raw haplotype data before filtering

**Quality Control:**
- `multiqc/multiqc_report.html` - Comprehensive quality metrics across all samples

**Sequoia Parentage Analysis (sequoia/ directory, if enabled):**
- `*_sequoia_parentage_results.txt` - Parent-offspring assignments
- `*_sequoia_output.rds` - Sequoia R object for troubleshooting or more detailed information.

**CKMRsim Analysis (if enabled):**
- `*_PO_results.tsv` - Parent-offspring assignments with quality metrics. Columns include:
  - `D2_indiv`: Offspring sample ID
  - `D1_indiv`: Parent sample ID
  - `logl_ratio`: Log-likelihood ratio for parent-offspring vs. unrelated
  - `num_loc`: Number of loci used in comparison
  - `total_incompat`: Number of Mendelian incompatibilities detected
  - `fraction_incompat`: Proportion of loci with incompatibilities (useful for quality control)
- `*_PO_ckmr.rds` - CKMRsim R object containing the CKMR model
- `*_offspring_genotypes_long.rds` - Processed offspring genotypes in long format
- `*_parent_genotypes_long.rds` - Processed parent genotypes in long format
- `*_summary.tsv` - Combined RUBIAS and CKMRsim results (when both are enabled). Found in summary subfolder within the main results folder.

## Troubleshooting

### Testing the Pipeline

A test configuration is available to verify your installation:

```bash
nextflow main.nf -profile docker -params-file test_settings.yml
```

The [test_settings.yml](test_settings.yml) file uses example data and can help diagnose configuration issues. For cloud testing, [test_settings_azure.yml](test_settings_azure.yml) provides an Azure-ready configuration. Fully documented template versions with parameter guidance can be found in [examples/settings/test_settings.yml](examples/settings/test_settings.yml) and [examples/settings/test_settings_azure.yml](examples/settings/test_settings_azure.yml). An equivalent JSON configuration is also available in [test_settings.json](test_settings.json).

### Common Issues

**Pipeline fails to start:**
- Verify Nextflow is installed: `nextflow -version`
- Check that Docker/Singularity is available, if using: `docker --version` or `singularity --version`
- Ensure your configuration files have correct paths for your environment
- For Azure profiles, ensure Azure storage credentials are configured in `~/.nextflow/config`

**Out of memory errors:**
- Increase memory allocation for specific processes in your profile config
- Check the `conf/` directory for examples of memory settings per process

**Missing reference files:**
- Verify all paths in your params file point to valid reference data
- For local/HPC profiles, ensure reference files are accessible from compute nodes

**Container pull failures:**
- Check network connectivity to Docker Hub or Singularity library
- For Singularity, verify the cache directory is writable (set via `NXF_SINGULARITY_CACHEDIR`)

### Getting Help

For additional support:
- Check the [GitHub Issues](https://github.com/DWR-CSI/chinook_gt/issues) page
- Review Nextflow logs in the `.nextflow.log` files
- Use `nextflow log` to view execution history

## Contributors

In addition to all contributors listed on GitHub, this code was originally adapted from a pipeline written by Jeff Rodzen and Joy Gaines at CDFW that utilizes code from the Garza Lab (NOAA SWFSC).

## License

This code is licensed for free non-commercial or commercial use under GPLv3.
