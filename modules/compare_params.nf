process COMPARE_PARAMS {
    tag "compare_params"

    conda './env.yaml'

    publishDir params.outdir, mode: 'copy'

    input:
    path(precision_recall_results)

    output:
    path "pr_curves.png" , emit: pr_curves
    path "barplot.png" , emit: barplot
    path "best_f1_scores.tsv" , emit: best_f1_scores
    stdout emit: out_string

    script:
    """
    compare_params.R ${precision_recall_results} best_f1_scores.tsv pr_curves.png barplot.png
    """
}
