process COMPARE_PARAMS {
    tag "compare_params"

    conda './env.yaml'

    publishDir "${params.top_outdir}/5_compare_params", mode: 'copy'

    input:
    path(precision_recall_results)

    output:
    path "barplot.png" , emit: barplot
    path "ADVANCED_pr_curves.png" , emit: pr_curves
    path "ADVANCED_scatterplot.png" , emit: scatterplot
    path "ADVANCED_best_f1_scores.tsv" , emit: best_f1_scores
    stdout emit: out_string

    script:
    """
    compare_params.R ${precision_recall_results} ADVANCED_best_f1_scores.tsv ADVANCED_pr_curves.png barplot.png ADVANCED_scatterplot.png
    """
}
