process COMPARE_PARAMS {
    tag "compare_params"

    publishDir params.outdir, mode: 'copy'

    input:
    path(precision_recall_results)

    output:
    path "pr_curves.png" , emit: pr_curves
    stdout emit: out_string

    script:
    """
    compare_params.R ${precision_recall_results} pr_curves.png
    """
}
