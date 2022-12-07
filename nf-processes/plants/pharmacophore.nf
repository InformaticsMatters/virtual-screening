params.interval = 1000
params.torsion_weight = 20.0
params.rmsd = 2.0
params.count = 10

process pharmacophore {

    container 'informaticsmatters/vs-plants:latest'

    input:
    path inputs // .sdf
    path fragments  // .mol or .sdf

    output:
    path "ph4_${inputs.name}"
    env COUNT

    """
    /code/pharmacophore.py\
      --input '$inputs'\
      --fragments '$fragments'\
      --outfile 'ph4_${inputs.name}'\
      --interval $params.interval\
      --count $params.count\
      --torsion-weight $params.torsion_weight\
      --rmsd $params.rmsd\
      --work-dir work

      # count the number of outputs
      COUNT=\$(fgrep -c '\$\$\$\$' 'ph4_${inputs.name}')
    """
}