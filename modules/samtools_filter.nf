process SAMTOOLS_FILTER {
    tag "$bam"

    conda './env.yaml'

    publishDir params.outdir, mode: 'copy'

    input:
    path bam

    output:
    path "${bam}.filtered", emit: filtered_bam

    script:
    """
    samtools view -b -F 2324 --threads ${task.cpus} ${bam} > ${bam}.filtered
    """
}
