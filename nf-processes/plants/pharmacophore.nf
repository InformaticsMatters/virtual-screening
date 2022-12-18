params.interval = 1000
params.torsion_weight = 20.0
params.rmsd = 2.0
params.threshold = 0
params.count = 10
params.gen3d = false

process pharmacophore {

    container 'informaticsmatters/vs-plants:latest'

    input:
    path inputs // .sdf or .smi
    path fragments  // .mol or .sdf

    output:
    path "ph4_${inputs.name}", optional: true
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
      ${params.threshold ? '--threshold ' + params.threshold : ''}\
      ${params.gen3d ? '--gen-coords' : ''}\
      --work-dir work

      # count the number of outputs - for some strange reason the fgrep command fails is the file is empty
      if [ -s 'o3da_${inputs.name}'
      then
        COUNT=\$(fgrep -c '\$\$\$\$' 'o3da_${inputs.name}')
      else
        COUNT=0
      fi
    """
}