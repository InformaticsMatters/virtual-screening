params.interval = 1000
params.output = null
params.score_mode = 'all'
params.tanimoto = false

process sucos {

    container 'informaticsmatters/vs-prep:3.0.0'

    input:
    path inputs // .sdf
    path reference  // .sdf or .mol

    output:
    path '*.sdf'
    path 'count.txt' // the number of outputs, used for the cost

    script:
    """
    OUT=${params.output == null ? 'sucos_' + inputs : params.output}
    /code/sucos.py\
      --input '$inputs'\
      --reference '$reference'\
      --output "\$OUT"\
      --interval $params.interval\
      --score-mode '$params.score_mode'\
      ${params.tanimoto ? '--tanimoto' : ''}\

      # count the number of outputs - for some strange reason the fgrep command fails is the file is empty
      if [ -s "\$OUT" ]
      then
        COUNT=\$(fgrep -c '\$\$\$\$' "\$OUT")
      else
        COUNT=0
      fi
      echo \$COUNT > count.txt
    """
}