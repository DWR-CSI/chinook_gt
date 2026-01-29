# chinook_gt

## Table of Contents
- [Overview](#overview)
- [Setup](#setup)
- [Usage](#usage)
- [Execution Profiles](#execution-profiles)
- [Parameters](#parameters)
  - [Input Options](#input-options)
  - [Trimmomatic Options](#trimmomatic-options)
  - [FLASH2 Options](#flash2-options)
  - [Processing Options](#processing-options)
  - [Microhaplotopia Options](#microhaplotopia-options)
  - [Threshold Options](#threshold-options)
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

Prior to running the primary `main.nf` workflow, simply ensure that the parameters and paths specified in the configuration files located at `nextflow.config`, `conf/azure.config` (or your preferred profile), and `~/.nextflow/config` are correct for your setup. You can then run `main.nf` with Nextflow while specifying a `.yml` parameters file corresponding to your run. See [examples/settings/test_settings.yml](examples/settings/test_settings.yml) for an example params-file.

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

- **`azure_local`** - Azure Batch with local adapter files
  - Similar to `azure` but uses local reference files instead of cloud storage

- **`docker`** - Local execution with Docker containers
  - For running on local machines or servers with Docker installed

- **`singularity`** - Singularity containers for HPC systems
  - Compatible with SLURM and other HPC schedulers
  - No root privileges required

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

### Input Options

- **`project`** (string, default: `"YourProject1"`)
  Project name (no spaces or special characters). Used as the output subdirectory name.

- **`input`** (path, default: `"az://seqera1/raw_data/ADPH_mushed/*_R{1,2}_concatenated.fastq.gz"`)
  Path to paired-end FastQ files.

### Trimmomatic Options

- **`trim_params`** (string, default: `"LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36"`)
  Trimmomatic adapter trimming parameters.

### FLASH2 Options

- **`min_overlap`** (integer, default: `3`)
  Minimum overlap length for read merging.

- **`max_overlap`** (integer, default: `400`)
  Maximum overlap length for read merging.

### Processing Options

- **`concat_all_reads`** (boolean, default: `false`)
  Concatenate all available reads per sample (merged, unmerged, and unpaired) for BWA mapping. If false, uses only merged and single-end reads.

- **`loci_to_remove`** (string, default: see below)
  Perl-compatible regular expression to match loci column names to remove from analysis. See [Loci Removal](#loci-removal) section for details.

### Full Genome Mapping Options

- **`genome_cache`** (path)
  Path to the genome cache directory. This directory is used to store the reference genome and BWA indices after they are downloaded and indexed for the first time. If not specified, the pipeline will attempt to use a default path of `${projectDir}/data/genome_cache`.

- **`full_genome_mount_path`** (path)
  **Required for Azure Batch execution (`-profile azure`)** when running the full panel in large numbers.
  Specifies the absolute path to the mounted Azure File Share directory containing the reference genome and BWA indices.
  Example: `/mnt/batch/tasks/fsmounts/seqres/reference/genomes`
  
  **Note:** The pipeline expects the genome file to be located at `Otsh_v1.0/Otsh_v1.0.fna` relative to this path, along with its associated BWA index files (e.g., `.bwt`, `.sa`, etc.).

### Microhaplotopia Options

- **`haplotype_depth`** (integer, default: `4`)
  Remove haplotypes with depth less than this value to filter potential genotyping errors.

- **`total_depth`** (integer, default: `8`)
  Minimum number of reads required for a genotype to be retained.

- **`allele_balance`** (number, default: `0.35`, range: 0-1)
  Minimum allele balance ratio required for genotypes to be retained. Calculated as the ratio of reads for a given microhaplotype divided by reads for the most common microhaplotype.

### Threshold Options

- **`ots28_missing_threshold`** (number, default: `0.7`, range: 0-1)
  If less than this proportion of OTS28 data is missing, consider OTS28 data Intermediate instead of uncertain.

- **`gsi_missing_threshold`** (number, default: `0.6`, range: 0-1)
  If more than this proportion of GSI data is missing, consider GSI data invalid.

- **`pofz_threshold`** (number, default: `0.8`, range: 0-1)
  If the maximum Probability of Z-score (PofZ) is less than this value, consider the result ambiguous.

### CKMRsim Options

  - **`use_CKMR`** (boolean, default: `false`)
  Enable CKMRsim for parent-offspring analysis. When enabled, the pipeline runs CKMRsim after genotype generation and can combine results with RUBIAS stock assignments.

  - **`CKMR_parent_geno_input`** (path, default: `"$projectDir/examples/PBT/FRH2024_reference_genotypes.csv"`)
  Path to parent reference genotype file in wide-format CSV. File format requirements:
  Header row with `SAMPLE_ID` as first column, followed by locus names with `.1` and `.2` suffixes (e.g., `Locus1.1`, `Locus1.2`) Missing data should be coded consistently (currently, we use NA)

  - **`CKMR_logl_threshold`** (number, default: `6.9`, minimum: 0)
  Log-likelihood ratio threshold for identifying parent-offspring relationships. Results with log-likelihood ratios below this threshold are filtered out.

  - **`CKMR_min_loci`** (integer, default: `90`, minimum: 1)
  Minimum number of loci required for a valid parent-offspring assignment.

  - **`CKMR_extra_genos_allele_freqs`** (path, default: `"$projectDir/examples/PBT/JPE2022-2024_geno_wide.csv"`)
  Path to extra genotypes file (wide-format CSV) containing additional individuals to include in allele frequency calculations. Including a larger and more representative sample improves the accuracy of kinship likelihood calculations. This file follows the same format as the parent genotype file.

### Sequoia Options

- **`use_sequoia`** (boolean, default: `false`)
  Enable Sequoia for parentage-based tagging (PBT) analysis.

- **`sequoia_mode`** (string, default: `"par"`, options: `"par"` or `"ped"`)
  Sequoia analysis mode: `"par"` for parentage analysis or `"ped"` for full pedigree analysis (slower).

- **`sequoia_missing_threshold`** (number, default: `0.5`, range: 0-1)
  If more than this proportion of data is missing, consider the result invalid. Set to 0.2 for conservative cutoff.

- **`parent_geno_input`** (path, default: `"az://seqera1/reference/full_panel/PBT/FRH2024_reference_genotypes.txt"`)
  Path to parent genotype file (Sequoia formatted).

- **`parent_lifehistory`** (path, default: `"az://seqera1/reference/full_panel/PBT/FRH2024_reference_lifehistory.txt"`)
  Path to parent metadata/lifehistory file (Sequoia formatted).

- **`offspring_birthyear`** (string, default: `"unknown"`)
  Birth year of the offspring.

- **`offspring_minBY`** (string, default: `"unknown"`)
  Minimum birth year for offspring.

- **`offspring_maxBY`** (string, default: `"2025"`)
  Maximum birth year for offspring. Set to current year if unknown.

- **`offspring_max_age`** (integer, default: `1`)
  Maximum age of the offspring in the dataset.

- **`species_max_repro_age`** (integer, default: `6`)
  Maximum reproductive age of the species.

- **`species_min_repro_age`** (integer, default: `1`)
  Minimum reproductive age of the species.

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

The [test_settings.yml](test_settings.yml) file uses example data and can help diagnose configuration issues.

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
