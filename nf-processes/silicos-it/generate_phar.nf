process generate_phar {

    container '3dechem/silicos-it:latest'

    input:
    file inputs

    output:
    file '*.phar'

    """
    align-it -d '$inputs' -p '${inputs.name[0..-5]}.phar'
    """
}
