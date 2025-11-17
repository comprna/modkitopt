#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modkitopt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline to optimise modkit pileup parameters --filter-threshold and --mod-threshold
    to maximise sensitivity and precision of nanopore RNA modification calls.
    
    Steps:
     1. Sort/filter modBAM with samtools
     2. Run modkit across parameter combinations
     3. Evaluate best parameters using validated (ground truth) sites
    
    Github : https://github.com/comprna/modkitopt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

nextflow.enable.dsl = 2

include { SAMTOOLS_SORT   } from './modules/samtools_sort.nf'
include { SAMTOOLS_FILTER } from './modules/samtools_filter.nf'
// include { MODKIT_PILEUP   } from './modules/modkit_pileup.nf'
// include { EVAL_PARAMS     } from './modules/eval_params.nf'

workflow {

    /*
     * =========================================================================
     *     Input validation
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

    // Stage files
    modbam       = file(params.modbam)
    ground_truth = file(ground_truth)
    
    /*
     * =========================================================================
     *     Sort and filter BAM
     * =========================================================================
     */
    
    // TODO: What happens if you don't do these steps?
    sorted_bam   = SAMTOOLS_SORT(modbam)
    filtered_bam = SAMTOOLS_FILTER(sorted_bam)
    // TODO: Do we need to index?
}
