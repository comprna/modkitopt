process MODKIT_PILEUP {
    tag "ft=${filter_threshold}, mt=${mod_threshold}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(bam), path(index), path(fasta), val(filter_threshold), val(mod_threshold)

    output:
    path "modkit_pileup_${filter_threshold}_${mod_threshold}.bed" , emit: modkit_bed
    path "modkit_pileup_${filter_threshold}_${mod_threshold}.log" , emit: modkit_log
    path "test_${filter_threshold}_${mod_threshold}.txt"

    script:
    """
    echo "${fasta} ${bam} ${filter_threshold} ${mod_threshold}" \\
    > "test_${filter_threshold}_${mod_threshold}.txt"

    ${params.modkit} pileup ${bam} \\
        "modkit_pileup_${filter_threshold}_${mod_threshold}.bed" \\
        --threads 8
        --filter-threshold A:${filter_threshold} \\
        --mod-threshold a:${mod_threshold} \\
        --log-filepath "modkit_pileup_${filter_threshold}_${mod_threshold}.log" \\
        --with-header \\
        --motif A 0 \\
        --ref ${fasta} \\
        --ignore 17596 # Remove inosine calls

    > "modkit_pileup_${filter_threshold}_${mod_threshold}.bed"
    """
}
