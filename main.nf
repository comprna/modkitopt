#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modkitopt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline to optimise the modkit pileup parameters --filter-threshold and 
    --mod-threshold to maximise the sensitivity and precision of nanopore RNA
    modification calls.
    
    Steps:
     1. Filter & sort modBAM with samtools
     2. Run modkit pileup across modkit parameter combinations
     3. For each modkit parameter combination, compute recall and precision of
        predicted sites using ground truth sites
     4. Determine modkit parameters that give highest recall and precision
    
    Github : https://github.com/comprna/modkitopt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

nextflow.enable.dsl = 2

include { SAMTOOLS_FILTER  } from './modules/samtools_filter.nf'
include { SAMTOOLS_SORT    } from './modules/samtools_sort.nf'
include { SAMTOOLS_INDEX   } from './modules/samtools_index.nf'
include { MODKIT_PILEUP    } from './modules/modkit_pileup.nf'
include { T2G              } from './modules/t2g.nf'
include { PRECISION_RECALL } from './modules/precision_recall.nf'
include { COMPARE_PARAMS   } from './modules/compare_params.nf'

def helpMessage() {
  log.info """
        Usage:
        The typical command structure for running the pipeline is as follows:
        nextflow run main.nf --modbam sample.bam
                             --mod_type <m6A|pseU|m5C|inosine>
                             --modkit /path/to/modkit
                             --fasta /path/to/transcriptome.fa
                             --annotation /path/to/annotation.gff3
                             -profile <local|pbs|pbspro|slurm>

        Mandatory arguments:
         --modbam             .bam file containing per-read modification calls
         --mod_type           Modification type (options: m6A, pseU, m5C, inosine)
         --modkit             Path to modkit executable
         --fasta              Path to reference transcriptome
         --annotation         Path to corresponding reference annotation (.gtf or .gff3)
         -profile             Execution environment (options: local, pbs, pbspro, slurm)

        Mandatory arguments if running on an HPC system (-profile is pbs, pbspro or slurm):
         --hpc_queue          Name of the queue that Nextflow can schedule jobs to (e.g., 'normal')
         --hpc_project        HPC project code that Nextflow can schedule jobs to (e.g., 'ab12')
         --hpc_storage        HPC storage location that outputs can be written to (e.g., 'gdata/ab12')
         --help               This usage statement

        Optional arguments:
         --truth_sites        .tsv file containing known modification sites (genomic coordinates, expected columns 1 and 2: [chr, pos], mandatory if mod_type is m5C or inosine)
        """
}

workflow {

    /*
     * =========================================================================
     *  Input validation
     * =========================================================================
     */
    
    if (params.help) {
      helpMessage()
      exit 0
    }
    if (!params.modbam) {
        exit 1, "ERROR: Please provide --modbam (.bam file output by modification caller)"
    }
    if (!params.mod_type) {
        exit 1, "ERROR: Please provide --mod_type (options: m6A, pseU, m5C, inosine)"
    }
    if (!params.modkit) {
        exit 1, "ERROR: Please provide the path to modkit via --modkit"
    }
    if (!params.fasta) {
        exit 1, "ERROR: Please provide the path to the reference transriptome via --fasta"
    }
    if (!params.annotation) {
        exit 1, "ERROR: Please provide the path to the reference annotation (.gtf or .gff3) via --annotation"
    }

    // HPC parameters must be specified if using pbs, pbspro or slurm
    def hpcProfiles = ['pbs','pbspro','slurm']

    if (workflow.profile in hpcProfiles) {
        if (!params.hpc_project) {
            exit 1, "ERROR: --hpc_project must be specified when using profile '${workflow.profile}'."
        }
        if (!params.hpc_storage) {
            exit 1, "ERROR: --hpc_storage must be specified when using profile '${workflow.profile}'."
        }
        if (!params.hpc_queue) {
            exit 1, "ERROR: --hpc_queue must be specified when using profile '${workflow.profile}'."
        }
    }

    // Validate modification type
    assert params.mod_type in ['m6A','pseU','m5C','inosine'], \
        "Modification type not recognised, choose from: m6A, pseU, m5C, inosine"

    // Validate ground truth
    def truth_sites
    if (!params.truth_sites) {

        // If no ground truth sites provided for m6A, then use supplied default
        if (params.mod_type == "m6A") {
            truth_sites = "./resources/m6A_validated.tsv"
            log.warn "No --truth_sites supplied for --mod_type m6A; using default '${truth_sites}'"

        // If no ground truth sites provided for pseU, then use supplied default
        } else if (params.mod_type == "pseU") {
            truth_sites = "./resources/pseU_validated.tsv"
            log.warn "No --truth_sites supplied for --mod_type pseU; using default '${truth_sites}'"

        // Ground truth sites must be provided for mods other than m6A or pseU
        } else {
            exit 1, "ERROR: --truth_sites is required when mod_type is not m6A or pseU"
        }

    } else {
        truth_sites = params.truth_sites
    }

    /*
     * =========================================================================
     *  Filter, sort and index BAM
     * =========================================================================
     */

    ch_modbam = channel.fromPath(params.modbam, checkIfExists: true)
    ch_filtered_bam = SAMTOOLS_FILTER(ch_modbam)
    ch_sorted_bam = SAMTOOLS_SORT(ch_filtered_bam)
    ch_indexed_bam = SAMTOOLS_INDEX(ch_sorted_bam)

    /*
     * =========================================================================
     *  Run modkit across parameter combinations
     * =========================================================================
     */

    // Create a channel of parameter combinations
    ch_modkit_params = channel.from(params.filter_thresholds)
                              .combine(channel.from(params.mod_thresholds))
                              .map { ft, mt -> tuple(ft, mt) }
                              .mix(Channel.value(["default", "default"]))

    // Create a channel for the reference fasta file
    ch_fasta = channel.fromPath(params.fasta, checkIfExists: true)

    // Combine the BAM and BAM index into a single channel
    ch_bam_index = ch_indexed_bam.indexed_bam.combine(ch_indexed_bam.index)

    // Combine also with the fasta channel
    ch_bam_fasta = ch_bam_index.combine(ch_fasta)

    // Combine also with the modkit parameters channel
    ch_modkit_input = ch_bam_fasta.combine(ch_modkit_params)

    // Run modkit for each parameter combination
    ch_modkit_output = MODKIT_PILEUP(ch_modkit_input)

    // Merge the modkit output bed file with the modkit parameters
    ch_bed_params = ch_modkit_output.modkit_bed
                                    .merge(ch_modkit_output.filter_threshold)
                                    .merge(ch_modkit_output.mod_threshold)

    /*
     * =========================================================================
     *  Convert transcriptomic to genomic coordinates for comparison with
     *  validated sites
     * =========================================================================
     */

    // Create a channel for the reference annotation
    ch_annotation = channel.fromPath(params.annotation, checkIfExists: true)

    // Combine the annotation with the bed file channel
    ch_t2g_input = ch_bed_params.combine(ch_annotation)

    // Convert transcriptomic to genomic coordinates
    ch_modkit_genomic = T2G(ch_t2g_input)

    // Merge the genomic bed file with the modkit parameters
    ch_genomic_params = ch_modkit_genomic.bed_genomic
                                         .merge(ch_modkit_genomic.filter_threshold)
                                         .merge(ch_modkit_genomic.mod_threshold)
    
    /*
     * =========================================================================
     *  Compute recall and precision of called sites, compared against ground
     *  truth sites, across stoichiometry thresholds
     * =========================================================================
     */
    
    // Create a channel for the ground truth sites
    ch_truth = channel.fromPath(truth_sites, checkIfExists: true)

    // Combine the channels for ground truth sites and modkit sites
    ch_eval_params_input = ch_genomic_params.combine(ch_truth)

    // Evaluate precision and recall
    ch_prec_recall = PRECISION_RECALL(ch_eval_params_input)

    /*
     * =========================================================================
     *  Compare performance of modkit parameters
     * =========================================================================
     */

    // Collect the results across all parameters
    ch_prec_recall_collected = ch_prec_recall.precision_recall.collect()

    ch_best_params = COMPARE_PARAMS(ch_prec_recall_collected)
    ch_best_params.out_string.view()

}
