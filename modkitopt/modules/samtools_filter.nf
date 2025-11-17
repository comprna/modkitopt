process SAMTOOLS_FILTER {
    tag "$sorted_bam"

    conda './env.yaml'

    input:
    path sorted_bam

    output:
    path "filtered.bam", emit: filtered_bam

    script:
    """
    samtools view -b -F 2324 ${sorted_bam} > filtered.bam
    """
}
