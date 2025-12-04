process SAMTOOLS_INDEX {
    tag "$bam"

    conda './env.yaml'

    publishDir params.outdir, mode: 'copy'

    input:
    path bam

    output:
    path "${bam}.bai" , emit: index
    path "${bam}" , emit: indexed_bam

    script:
    """
    samtools index ${bam}
    """
}
