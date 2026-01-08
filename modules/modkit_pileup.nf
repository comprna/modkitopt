process MODKIT_PILEUP {
    tag "ft=${filter_threshold}, mt=${mod_threshold}"

    publishDir "${params.top_outdir}/2_pileup", mode: 'copy'

    input:
    tuple path(bam),
          path(index),
          path(fasta),
          val(filter_threshold),
          val(mod_threshold),
          val(mod_type)

    output:
    path "modkit_pileup_${filter_threshold}_${mod_threshold}.bed" , emit: modkit_bed
    path "modkit_pileup_${filter_threshold}_${mod_threshold}.log" , emit: modkit_log
    val(filter_threshold) , emit: filter_threshold
    val(mod_threshold) , emit: mod_threshold

    script:

    def canonical_base_map = [
      m6A:     'A',
      pseU:    'T',
      m5C:     'C',
      inosine: 'A'
    ]

    def mod_code_map = [
      m6A:     'a',
      pseU:    '17802',
      m5C:     'm',
      inosine: '17596'
    ]

    def canonical_base = canonical_base_map[mod_type]
    def mod_code = mod_code_map[mod_type]

    def filter_opt = filter_threshold == 'default' ? '' :
                      "--filter-threshold ${canonical_base}:${filter_threshold}"

    def mod_opt    = mod_threshold == 'default' ? '' :
                      "--mod-threshold ${mod_code}:${mod_threshold}"

    """
    ${params.modkit} pileup ${bam} \\
        "modkit_pileup_${filter_threshold}_${mod_threshold}.bed" \\
        --threads ${task.cpus} \\
        --modified-bases ${mod_type} \\
        --with-header \\
        ${filter_opt} \\
        ${mod_opt} \\
        --log-filepath "modkit_pileup_${filter_threshold}_${mod_threshold}.log" \\
        --ref ${fasta} \\
        --preload-references
    """
}
