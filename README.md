# ModkitOpt

ModkitOpt finds the best `--mod-thresold` and `--filter-threshold` parameters to use when running `modkit pileup`, and the best stoichiometry for filtering modkit's bedMethyl output, to maximise the sensitivity and recall of your nanopore direct RNA modification calls.

### Why use ModkitOpt?

By default `--mod-threshold` and `--filter-threshold` are estimated by modkit using an algorithm that does not account for ... , which can result in false modification calls.

### How ModkitOpt works

ModkitOpt ...

Note on providing ground truth sites


### Citation

If you use this software, please cite:

> Ref TBD


# Quick start

To run locally using the example modBAM file we provide in modkitopt/resources, simply:

**1. Clone the repository**

```
git clone ...
```

**2. Install dependencies**

* modkit >= v0.6.0

  * Download modkit_vXYZ.tar.gz from the [modkit release page](https://github.com/nanoporetech/modkit/releases)
  * Extract the archive contents
  ```
  tar -xvzf modkit_vXYZ.tar.gz
  ```

* nextflow
* conda

Nextflow will automatically install all other dependencies using conda (environment defined in `modkitopt/env.yaml`) the first time that modkitopt is run.

**3. Run modkitopt**

```bash
cd /path/to/modkitopt

nextflow run main.nf \
  --modbam ./resources/example.bam \
  --mod_type m6A \
  --modkit /path/to/modkit \
  --fasta ./resources/gencode.v45.transcripts.fa \
  --gff ./resources/gencode.v45.annotation.gff3 \
  --ground_truth ./resources/m6A_validated.tsv \
  -profile local
```

*Note:* The first time that you run modkitopt, nextflow will create a conda environment and install dependencies - be patient, this will take a few minutes.

# Running in HPC environments

Before running modkitopt inside a job, first run it on a login node (or a node where internet is available) so that nextflow can create the conda environment.

## Dependency installation

Most of the following dependencies are often already provided in HPC environments as modules that can simply be loaded. If not, they need to be installed following the guidelines for your system.

* modkit >= v0.6.0 (copy the downloaded binary to your HPC storage location)
* nextflow >= X
* samtools
* python3
* R >= X

### R packages

```bash
$ R

> install.packages(c("glue", "tidyverse", "rlist", "ggrepel","BiocManager"))
> BiocManager::install(c("GenomicFeatures", "txdbmaker"))
```

*Note:* To enable R package installation in your HPC environment, you may need to first install some intel compiler packages. See your HPC documentation for guidance (e.g. For NCI's gadi: https://opus.nci.org.au/spaces/Help/pages/248840714/R).

### A note on containers and conda in HPC environments

We don't provide a container for running modkitopt since containers are not ideal for HPC environments. Although a container can be run anywhere, the installed software binaries are built to run on the container’s operating system and are therefore not optimised for the HPC environment's architecture, therefore potentially compromising performance. Similarly, conda is disabled by default for HPC environments (although you can modify this in nextflow.config by setting conda.enabled = true in the relevant profile).


## Run ModkitOpt

To enable nextflow to ... specify your HPC environment with the `-profile` flag, such as `-profile pbs` or `-profile slurm`. For NCI's gadi use `-profile pbspro`


```
Command to run
```

Overriding parameters:
nextflow run main.nf -profile pbspro \
    --pileup_cpus 4 \
    --pileup_memory 16GB
    --hpc_project
    --hpc_storage


### Resuming an interrupted run

If your run gets interrupted, Nextflow automatically supports checkpointing and resuming runs. Simply add `-resume` to the Nextflow command that didn't complete and run again!

For example:

```bash
nextflow run main.nf \
  --modbam ./resources/example.bam \
  --mod_type m6A \
  --modkit /path/to/modkit \
  --fasta ./resources/gencode.v45.transcripts.fa \
  --gff ./resources/gencode.v45.annotation.gff3 \
  --ground_truth ./resources/m6A_validated.tsv \
  -profile local
  -resume
```

### How resources are allocated

-- Explain Nextflow manages threads (for pbs, slurm then 8 threads allocated per modkit pileup call - maximum recommended by ONT)

### Estimated run-time

-- Estimated run-time (using examples: HPC environment, local computer, human transcriptome ref)

# Tested environments and software versions

-- What versions of everything has been tested (nextflow version, modkit version...)