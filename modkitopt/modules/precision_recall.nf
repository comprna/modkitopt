process PRECISION_RECALL {
    tag "ft=${filter_threshold}, mt=${mod_threshold}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(bed_genomic),
          val(filter_threshold),
          val(mod_threshold),
          path(ground_truth)

    output:
    path "out_${filter_threshold}_${mod_threshold}.txt"

    script:
    """
    echo ${bed_genomic} ${ground_truth} ${filter_threshold} ${mod_threshold} \
        > out_${filter_threshold}_${mod_threshold}.txt
    """
}



    // python3 evaluate_params.py \\
    //     --truth ${truth} \\
    //     --results modkit_results.pkl \\
    //     --out-best best_parameters.txt \\
    //     --out-summary evaluation_summary.tsv

    // python3 eval_recall_precision.py \\
    //     --input_bed ${bed_genomic}