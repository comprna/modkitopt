process MODKIT_PILEUP {
    tag "ft=${filter_threshold}, mt=${mod_threshold}"

    publishDir params.outdir, mode: 'copy'

    cpus params.pileup_cpus ?: 1

    input:
    tuple path(bam),
          path(index),
          path(fasta),
          val(filter_threshold),
          val(mod_threshold)

    output:
    path "modkit_pileup_${filter_threshold}_${mod_threshold}.bed" , emit: modkit_bed
    path "modkit_pileup_${filter_threshold}_${mod_threshold}.log" , emit: modkit_log
    val(filter_threshold) , emit: filter_threshold
    val(mod_threshold) , emit: mod_threshold

    script:
    """
    ${params.modkit} pileup ${bam} \\
        "modkit_pileup_${filter_threshold}_${mod_threshold}.bed" \\
        --modified-bases ${params.mod_type} \\
        --threads ${task.cpus} \\
        --with-header \\
        --filter-threshold A:${filter_threshold} \\
        --mod-threshold a:${mod_threshold} \\
        --log-filepath "modkit_pileup_${filter_threshold}_${mod_threshold}.log" \\
        --ref ${fasta} \\
    """
}
