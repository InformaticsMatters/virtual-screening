params.interval = 1000
params.crippen = false
params.remove_hydrogens = true

process open3dalign {

    container 'informaticsmatters/vs-prep:3.1.0'

    input:
    path inputs // .sdf
    path query  // .sdf or .mol
    val threshold // 0 or null for no threshold

    output:
    path "o3da_${inputs.name}"
    path 'count.txt' // the number of outputs, used for the cost

    script:
    """
    /code/open3dalign.py\
      --inputs '$inputs'\
      --query '$query'\
      --outfile 'o3da_${inputs.name}'\
      --interval $params.interval\
      ${params.remove_hydrogens ? '--remove-hydrogens' : ''}\
      ${params.crippen ? '--crippen' : ''}\
      ${threshold ? '--threshold ' + threshold : ''}

      # count the number of outputs - for some strange reason the fgrep command fails is the file is empty
      if [ -s 'o3da_${inputs.name}' ]
      then
        COUNT=\$(fgrep -c '\$\$\$\$' 'o3da_${inputs.name}')
      else
        COUNT=0
      fi
      echo \$COUNT > count.txt
    """
}