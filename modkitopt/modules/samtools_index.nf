process SAMTOOLS_INDEX {
    tag "$bam"

    conda './env.yaml'

    publishDir 'results', mode: 'copy'

    input:
    path bam

    output:
    path "${bam}.bai"
    path "${bam}" , emit: indexed_bam

    script:
    """
    samtools index ${bam}
    """
}
