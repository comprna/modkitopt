process T2G {
    tag "${bed}"

    conda './env.yaml'

    publishDir params.outdir_t2g, mode: 'copy'

    input:
    tuple path(bed), val(filter_threshold), val(mod_threshold), path(gff)

    output:
    path "${bed}.genomic" , emit: bed_genomic
    val(filter_threshold) , emit: filter_threshold
    val(mod_threshold) , emit: mod_threshold

    script:
    """
    t2g.R ${bed} ${gff}
    """
}
