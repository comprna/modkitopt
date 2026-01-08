process INSTALL_ENV {
    tag "install"

    conda './env.yaml'

    input:
    val(dummy)

    output:
    val(true)

    script:
    """
    echo "Conda environment successfully installed"
    """
}
