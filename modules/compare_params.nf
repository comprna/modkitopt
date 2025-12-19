process COMPARE_PARAMS {
    tag "compare_params"

    conda './env.yaml'

    publishDir "${params.top_outdir}/5_compare_params", mode: 'copy'

    input:
    path(precision_recall_results)

    output:
    path "barplot.png" , emit: barplot
    path "advanced/pr_curves.png" , emit: pr_curves
    path "advanced/scatterplot.png" , emit: scatterplot
    path "advanced/best_f1_scores.tsv" , emit: best_f1_scores
    stdout emit: out_string

    script:
    """
    compare_params.R ${precision_recall_results} advanced/best_f1_scores.tsv advanced/pr_curves.png barplot.png advanced/scatterplot.png
    """
}
