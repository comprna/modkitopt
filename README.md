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

* modkit >= v0.6.0

  * Download modkit_vXYZ.tar.gz from the [modkit release page](https://github.com/nanoporetech/modkit/releases)
  * Extract the archive contents
    ```
    tar -xvzf modkit_vXYZ.tar.gz
    ```
* conda ([Miniconda installation guide](https://www.anaconda.com/docs/getting-started/miniconda/install))
* nextflow ([installation guide](https://www.nextflow.io/docs/latest/install.html))


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

## Dependencies

The **modkit** binary that you downloaded and extracted in [Quick start](#quick-start) can simply be copied to your HPC storage location.

**Nextflow** and **conda** are often already provided in HPC environments as modules that can simply be loaded. If not, they need to be installed following the guidelines for your system.

## Run ModkitOpt

Before running modkitopt inside a job, first run it on a login node (or a node where internet is available) so that nextflow can create the conda environment. Once the conda environment is created and the nextflow pipeline starts executing, you can kill the pipeline and then start it in a job.

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


### Estimated run-time

-- Estimated run-time (using examples: HPC environment, local computer, human transcriptome ref)

# Tested environments and software versions

-- What versions of everything has been tested (nextflow version, modkit version...)