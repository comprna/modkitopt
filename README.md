# ModkitOpt

ModkitOpt finds the best `--mod-thresold` and `--filter-threshold` parameters to use when running `modkit pileup`, and the best stoichiometry for filtering modkit's bedMethyl output, to maximise the sensitivity and recall of your nanopore direct RNA modification calls.

### Why use ModkitOpt?

By default, modkit estimates `--mod-threshold` and `--filter-threshold` by considering the confidence of dorado per-read modification predictions, without any knowledge of whether prediction confidence correlates with prediction accuracy. This means that for datasets where modification predictions overall have very low confidence (as is often the case for rare modification types), the modkit thresholds in turn will be low, which allows low confidence calls to pass the threshold and therefore lead to sites being incorrectly predicted as modified. Conversely, high thresholds can result in missing true sites. We show in our paper (referenced below) that the thresholds estimated by modkit are not optimised for precision or sensitivity.

### How ModkitOpt works

ModkitOpt ...

Note on providing ground truth sites


### Citation

If you use this software, please cite:

> Ref TBD


# Quick start

To run locally using the example modBAM file we provide in `modkitopt/resources`, simply:

**1. Clone the repository**

```
git clone ...
```

**2. Install dependencies**

* **modkit** >= v0.6.0

  * Download modkit_vXYZ.tar.gz from the [modkit release page](https://github.com/nanoporetech/modkit/releases)
  * Extract the archive contents
    ```
    tar -xvzf modkit_vXYZ.tar.gz
    ```
* **conda** ([Miniconda installation guide](https://www.anaconda.com/docs/getting-started/miniconda/install))
* **nextflow** ([installation guide](https://www.nextflow.io/docs/latest/install.html))


Nextflow will automatically install all other dependencies using conda (environment defined in `modkitopt/env.yaml`) the first time that modkitopt is run.

**3. Run ModkitOpt**

```bash
cd /path/to/modkitopt

nextflow run main.nf                                          \
  --modbam           ./resources/example.bam                  \
  --mod_type         m6A                                      \
  --modkit           /path/to/modkit                          \
  --fasta            ./resources/gencode.v45.transcripts.fa   \
  --gff              ./resources/gencode.v45.annotation.gff3  \
  --ground_truth     ./resources/m6A_validated.tsv            \
  -profile local
```

*Note:* The first time that you run ModkitOpt, Nextflow will create a conda environment and install dependencies - be patient, this will take a few minutes.

# Running in HPC environments

We recommend running ModkitOpt in an HPC environment, since modkit is called multiple times with different threshold values. Nextflow handles submitting modkit jobs so that they can run at the same time to reduce the overall execution time of ModkitOpt.

## Dependencies

The **modkit** binary that you downloaded and extracted in [Quick start](#quick-start) can simply be copied to your HPC storage location.

**Nextflow** and **conda** are often already provided in HPC environments as modules that can simply be loaded. If not, they need to be installed following the guidelines for your system.

## Running ModkitOpt

### The first time you run ModkitOpt in an HPC environment

Before running ModkitOpt inside a job, first run it on a login node (or a node where internet is available) so that nextflow can create the conda environment. Once the conda environment is created and the nextflow pipeline starts executing, you can kill the pipeline and then proceed with submitting your Nextflow job.

### Specifying your HPC environment details

When running in an HPC environment, you need to specify three things:

**1. Your HPC environment profile**

This tells Nextflow what type of workload manager it is dealing with. We currently support PBS, PBS Pro and Slurm systems. Specify this with the `-profile` flag, such as `-profile pbs`, `-profile pbspro` or `-profile slurm`. For NCI's gadi use `-profile pbspro`. Nextflow automatically handles creating and submitting jobs in each of these environments.

*Note:* We have only tested ModkitOpt in a `pbspro` environment (NCI's gadi). While we have written profiles for `pbs` and `slurm`, these have not been tested. We welcome contributions from the community to improve these profiles, which can be found in `nextflow.config`.

**2. Your HPC project code**

You must specify the HPC project code that Nextflow can schedule jobs to using the `--hpc_project` flag, such as `--hpc_project ab12`.

**3. Your HPC storage location**

You must specify your HPC storage location using the `--hpc_storage` flag, such as `--hpc_storage gdata/ab12`.

### Command example

Briefly, the required input files are:
1. ModBAM file output by dorado
2. FASTA and GFF file for modkit to use (the same FASTA as you provided to dorado, with corresponding GFF annotation).
3. TSV file containing ground truth sites (optional if you're working with human data and your modification type is m6A or pseU)

See [Command details](#command-details) for more information.

```bash
nextflow run main.nf                                           \
  --modbam          /path/to/modbam.bam                        \
  --mod_type        m6A                                        \
  --modkit          /path/to/dist_modkit_v0.6.0_68b540b/modkit \
  --fasta           /path/to/ref.fa                            \
  --gff             /path/to/annotation.gff3                   \
  --ground_truth    /path/to/ground_truth_sites.tsv            \
  -profile          pbspro                                     \
  --hpc_project     ab12                                       \
  --hpc_storage     gdata/ab12+gdata/cd34
```


### Resuming an interrupted run

If your run gets interrupted, Nextflow automatically supports checkpointing and resuming runs. Simply add `-resume` to the Nextflow command that didn't complete and run again!

For example:

```bash
nextflow run main.nf                                           \
  --modbam          /path/to/modbam.bam                        \
  --mod_type        m6A                                        \
  --modkit          /path/to/dist_modkit_v0.6.0_68b540b/modkit \
  --fasta           /path/to/ref.fa                            \
  --gff             /path/to/annotation.gff3                   \
  --ground_truth    /path/to/ground_truth_sites.tsv            \
  -profile          pbspro                                     \
  --hpc_project     ab12                                       \
  --hpc_storage     gdata/ab12+gdata/cd34                      \
  -resume
```


### Estimated run-time

-- Estimated run-time (using examples: HPC environment, local computer, human transcriptome ref)


### Advanced use: Overriding the default resource allocation



    --pileup_cpus 4 \
    --pileup_memory 16GB
    default_cpus      = 1               // Default CPU count to request per process
    default_memory    = '8GB'           // Default memory to request per process
    default_time      = '1h'            // Default time to request per process
    pileup_cpus       = 8               // Default CPU count to request for modkit pileup
    pileup_memory     = '8GB'           // Default memory to request for modkit pileup
    pileup_time       = '4h'            // Default time to request for modkit pileup


# Command details



# Tested environments and software versions

-- What versions of everything has been tested (nextflow version, modkit version...)