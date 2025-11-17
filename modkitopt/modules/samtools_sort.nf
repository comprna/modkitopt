process SAMTOOLS_SORT {
    tag "$modbam"

    conda './env.yaml'

    input:
    path modbam

    output:
    path "sorted.bam" , emit: sorted_bam

    script:
    """
    samtools sort ${modbam} > "sorted.bam"
    echo "samtools sort: input ${modbam}, emitting sorted.bam"
    """
}
