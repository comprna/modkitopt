process MODKIT_PILEUP {
    tag "ft=${filter_threshold}, mt=${mod_threshold}"

    publishDir params.outdir_pileup, mode: 'copy'

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

    def filter_opt = filter_threshold == 'default' ? '' :
                      "--filter-threshold A:${filter_threshold}"

    def mod_opt    = mod_threshold == 'default' ? '' :
                      "--mod-threshold a:${mod_threshold}"

    """
    ${params.modkit} pileup ${bam} \\
        "modkit_pileup_${filter_threshold}_${mod_threshold}.bed" \\
        --threads ${task.cpus} \\
        --modified-bases ${params.mod_type} \\
        --with-header \\
        ${filter_opt} \\
        ${mod_opt} \\
        --log-filepath "modkit_pileup_${filter_threshold}_${mod_threshold}.log" \\
        --ref ${fasta} \\
        --preload-references
    """
}
