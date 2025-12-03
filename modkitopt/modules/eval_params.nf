process EVAL_PARAMS {
    tag "param_eval"

    input:
    tuple path(bed_genomic), path(params.ground_truth)

    output:

    script:
    """
    echo ${bed_genomic}
    echo ${params.ground_truth}
    """
}



    // python3 evaluate_params.py \\
    //     --truth ${truth} \\
    //     --results modkit_results.pkl \\
    //     --out-best best_parameters.txt \\
    //     --out-summary evaluation_summary.tsv

    // python3 eval_recall_precision.py \\
    //     --input_bed ${bed_genomic}