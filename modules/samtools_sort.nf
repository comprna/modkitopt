process SAMTOOLS_SORT {
    tag "$bam"

    conda './env.yaml'

    publishDir "${params.top_outdir}/samtools", mode: 'copy'

    input:
    path bam

    output:
    path "${bam}.sorted" , emit: sorted_bam

    script:
    """
    samtools sort --threads ${task.cpus} ${bam} > ${bam}.sorted
    """
}
