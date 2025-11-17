process SAMTOOLS_FILTER {
    tag "$sorted_bam"

    input:
    path sorted_bam

    output:
    path "filtered.bam", emit: filtered_bam

    conda "samtools"
    container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    script:
    """
    samtools view -F 0x900 -b ${sorted_bam} > filtered.bam
    """
}
