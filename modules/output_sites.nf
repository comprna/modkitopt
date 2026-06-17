process OUTPUT_SITES {
    tag "output_sites"

    conda './env.yaml'

    publishDir "${params.top_outdir}/5_outputs", mode: 'copy'

    input:
    path(best_params)
    path(modkit_beds)

    output:
    path "optimal_thresholds_modkit_pileup_*.bed" , emit: methylbed

    script:
    """
    filter_threshold=\$(awk 'NR==2 {print \$1}' ${best_params})
    mod_threshold=\$(awk 'NR==2 {print \$2}' ${best_params})

    methylbed="modkit_pileup_\${filter_threshold}_\${mod_threshold}.bed"
    selected="optimal_thresholds_\${methylbed}"

    echo "Copying: \${methylbed}"

    cp "\${methylbed}" "\${selected}"
    """
}
