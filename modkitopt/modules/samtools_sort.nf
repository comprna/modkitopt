process SAMTOOLS_SORT {
    tag "$bam"

    conda './env.yaml'

    publishDir 'results', mode: 'copy'

    input:
    path bam

    output:
    path "${bam}.sorted" , emit: sorted_bam

    script:
    """
    samtools sort ${bam} > ${bam}.sorted
    echo "samtools sort: input ${bam}, emitting ${bam}.sorted"
    """
}
