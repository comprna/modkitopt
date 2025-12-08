process PRECISION_RECALL {
    tag "ft=${filter_threshold}, mt=${mod_threshold}"

    conda './env.yaml'

    publishDir params.outdir_pr, mode: 'copy'

    input:
    tuple path(bed_genomic),
          val(filter_threshold),
          val(mod_threshold),
          path(ground_truth)

    output:
    path "precision_recall_${filter_threshold}_${mod_threshold}.tsv" , emit: precision_recall
    val(filter_threshold) , emit: filter_threshold
    val(mod_threshold) , emit: mod_threshold

    script:
    """
    precision_recall.py \
        --input_bed ${bed_genomic} \
        --truth ${ground_truth} \
        --output precision_recall_${filter_threshold}_${mod_threshold}.tsv
    """
}
