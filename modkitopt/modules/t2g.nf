process T2G {
    tag "${bed}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(bed), path(gff)

    output:
    path "${bed}.genomic" , emit: bed_genomic

    script:
    """
    Rscript ${workflow.projectDir}/bin/t2g.R ${bed} ${gff}
    """
}
