# chinook_gt

## Overview

This Nextflow pipeline takes .FASTQ files as input and produces Chinook run identifications as output along with relevant intermediate files. 

## Setup
Currently, this pipeline is optimized for use with Microsoft Azure Batch (and other clouds) but can/will be adapted for use with HPC systems by specifying Singularity containers instead of Docker containers.

Prior to running the primary `main.nf` workflow, simply ensure that the parameters and paths specified in the configuration files located at `nextflow.config`, `conf/azure.config` (or your preferred profile), and `~/.nextflow/config` are correct for your setup. You can then run `main.nf` with Nextflow while specifying a `.yml` parameters file corresponding to your run. See [examples/settings/full_settings.yml](settings.yml) for an example params-file.

### Loci Removal

The pipeline supports removing specific loci from analysis using the `loci_to_remove` parameter. This uses Perl-compatible regular expressions to match loci column names.

**How it works:**
- The regex pattern is matched against these column names in both RUBIAS and Sequoia analyses
- Matched loci are removed before Rubias and Sequoia analyses begins, but after the point of haplotype/genotype generation.

**Configuration:**
In your parameters YAML file:
```yaml
loci_to_remove: "(NC_037130\\.1:1062935-1063235)|(NC_037130\\.1:864908-865208)|(NC_037104\\.1:56552952-56553042)"
```

**Notes:**
- Leave empty (`""`) to remove no loci
- Use Perl-compatible regex syntax
- Cannot contain single quotes
- Backslashes must be escaped (use `\\.` to match a literal period)
- Invalid patterns will be caught and no loci will be removed

Example usage:

```
nextflow main.nf -profile azure -params-file settings.yml
```

If running locally instead of on Azure:
```
nextflow main.nf -profile docker -params-file settings.yml
```

## Output
The workflow will produce various intermediate and final files in the output directory that you specify. Final genetic run IDs can be found in the rubias subfolder. The file ending in `_summary.tsv` should be used for DWR Central Valley Chinook salmon run identification because it includes an RoSA-informed genetic run ID.

## Contributors

In addition to all contributors listed on GitHub, this code was originally adapted from a pipeline written by Jeff Rodzen and Joy Gaines at CDFW that utilizes code from the Garza Lab (NOAA SWFSC).

## License

This code is licensed for free non-commercial or commercial use under GPLv3.