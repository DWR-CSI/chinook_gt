# chinook_gt

## Overview

This Nextflow pipeline takes .FASTQ files as input and produces Chinook run identifications as output along with relevant intermediate files. 

## Setup
Currently, this pipeline is optimized for use with Microsoft Azure Batch (and other clouds) but can/will be adapted for use with HPC systems by specifying Singularity containers instead of Docker containers.

Prior to running the primary `main.nf` workflow, simply ensure that the parameters and paths specified in the configuration files located at `nextflow.config`, `conf/azure.config` (or your preferred profile), and `~/.nextflow/config` are correct for your setup. You can then run `main.nf` with Nextflow while specifying a `.yml` parameters file corresponding to your run. See [settings.yml](settings.yml) for an example params-file.

Example usage:

```
nextflow main.nf -profile azure -params-file settings.yml
```

## Contributors

In addition to all contributors listed on GitHub, this code was originally adapted from a pipeline written by Jeff Rodzen and Joy Gaines at CDFW that utilizes code from the Garza Lab (NOAA SWFSC).

## License

This code is licensed for free non-commercial or commercial use under GPLv3.