process EVAL_PARAMS {
    tag "param_eval"

    input:
    path truth from params.truth
    tuple val(filter), val(modthresh), path(out) from modkit_results.collect()

    output:
    path "best_parameters.txt"
    path "evaluation_summary.tsv"

    script:
    """
    python3 evaluate_params.py \\
        --truth ${truth} \\
        --results modkit_results.pkl \\
        --out-best best_parameters.txt \\
        --out-summary evaluation_summary.tsv
    """
}