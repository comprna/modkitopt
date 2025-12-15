process SAMTOOLS_INDEX {
    tag "$bam"

    conda './env.yaml'

    publishDir "${params.top_outdir}/samtools", mode: 'copy'

    input:
    path bam

    output:
    path "${bam}.bai" , emit: index
    path "${bam}" , emit: indexed_bam

    script:
    """
    samtools index --threads ${task.cpus} ${bam}
    """
}
