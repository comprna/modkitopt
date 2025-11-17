process SAMTOOLS_SORT {
    // Label to identify process being executed
    tag "$modbam"

    input:
    path modbam from params.modbam

    output:
    path "sorted.bam", emit: sorted_bam

    script:
    """
    #samtools sort -o sorted.bam ${modbam}
    echo "Samtools sort, input ${modbam}, emitting ${sorted_bam}"
    """
}
