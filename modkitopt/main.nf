#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modkitopt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline to optimise the modkit pileup parameters --filter-threshold and 
    --mod-threshold to maximise the sensitivity and precision of nanopore RNA
    modification calls.
    
    Steps:
     1. Sort/filter modBAM with samtools
     2. Run modkit across parameter combinations
     3. Evaluate best parameters using validated (ground truth) sites
    
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

workflow {

    /*
     * =========================================================================
     *  Input validation
     * =========================================================================
     */
    
    if (!params.modbam) {
        exit 1, "ERROR: Please provide --modbam (.bam file output by modification caller)"
    }
    if (!params.mod_type) {
        exit 1, "ERROR: Please provide --mod_type (options: m6A, m5C, pseU, other)"
    }
    if (!params.modkit) {
        exit 1, "ERROR: Please provide the path to modkit via --modkit"
    }
    if (!params.fasta) {
        exit 1, "ERROR: Please provide the path to the reference transriptome via --fasta"
    }
    if (!params.gff) {
        exit 1, "ERROR: Please provide the path to the reference annotation via --gff"
    }


    // Validate modification type
    assert params.mod_type in ['m6A','m5C','pseU','other'], \
        "Modification type not recognised, choose from: m6A, m5C, pseU, other"

    // Validate ground truth
    def ground_truth
    if (!params.ground_truth) {

        // If no ground_truth provided for m6A, then use supplied default
        if (params.mod_type == "m6A") {
            ground_truth = "./resources/m6A_validated.pickle"
            log.warn "No --ground_truth supplied for --mod_type m6A; using default '${ground_truth}'"

        // ground_truth must be provided for mods other than m6A
        } else {
            exit 1, "ERROR: --ground_truth is required when mod_type is not m6A"
        }

    } else {
        ground_truth = params.ground_truth
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
    ch_gff = channel.fromPath(params.gff, checkIfExists: true)

    // Combine the annotation with the bed file channel
    ch_t2g_input = ch_bed_params.combine(ch_gff)

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
    ch_truth = channel.fromPath(ground_truth, checkIfExists: true)

    // Combine the channels for ground truth sites and modkit sites
    ch_eval_params_input = ch_genomic_params.combine(ch_truth)

    // Evaluate precision and recall
    ch_prec_recall = PRECISION_RECALL(ch_eval_params_input)

}


// /home/alex/Documents/tools/dist_modkit_v0.5.1_8fa79e3/modkit
// /mnt/sda/projects/m6A_proteins/1_prelim_analysis/moved_from_OneDrive/1_prelim_analysis/0_refs/gencode/gencode.v45.transcripts.fa