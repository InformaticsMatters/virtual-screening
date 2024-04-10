process generate_phar {

    // latest is the only tag they have. Not been updated for 7 years.
    container '3dechem/silicos-it:latest'

    input:
    file inputs

    output:
    file '*.phar'

    """
    align-it -d '$inputs' -p '${inputs.name[0..-5]}.phar'
    """
}
