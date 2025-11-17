process MODKIT_PILEUP {
    tag "ft=${filter_threshold}, mt=${mod_threshold}"

    input:
    path bamfile from cleaned_bam
    val modtype from validated_modtype
    tuple val(filter_threshold), val(mod_threshold) from modkit_param_grid

    output:
    tuple val(filter_threshold), val(mod_threshold), path("modkit_${filter_threshold}_${mod_threshold}.tsv") \
        into modkit_results

    script:
    """
    modkit pileup \\
        --filter-threshold ${filter_threshold} \\
        --mod-threshold ${mod_threshold} \\
        --only-mod ${modtype} \\
        ${bamfile} \\
        modkit_${filter_threshold}_${mod_threshold}.tsv
    """
}