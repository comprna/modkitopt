process PRECISION_RECALL {
    tag "ft=${filter_threshold}, mt=${mod_threshold}"

    conda './env.yaml'

    publishDir "${params.top_outdir}/4_precision_recall", mode: 'copy'


    input:
    tuple path(bed_genomic),
          val(filter_threshold),
          val(mod_threshold),
          path(truth_sites)

    output:
    path "precision_recall_${filter_threshold}_${mod_threshold}.tsv" , emit: precision_recall
    val(filter_threshold) , emit: filter_threshold
    val(mod_threshold) , emit: mod_threshold

    script:
    """
    precision_recall.py \
        --input_bed ${bed_genomic} \
        --truth ${truth_sites} \
        --output precision_recall_${filter_threshold}_${mod_threshold}.tsv
    """
}
