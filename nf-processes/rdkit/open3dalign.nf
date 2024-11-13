params.interval = 1000
params.crippen = false
params.remove_hydrogens = true
params.threshold = 0

process open3dalign {

    container 'informaticsmatters/vs-prep:latest'

    input:
    path inputs // .sdf
    path query  // .sdf or .mol

    output:
    path "o3da_${inputs.name}"
    env COUNT

    """
    /code/open3dalign.py\
      --inputs '$inputs'\
      --query '$query'\
      --outfile 'o3da_${inputs.name}'\
      --interval $params.interval\
      ${params.remove_hydrogens ? '--remove-hydrogens' : ''}\
      ${params.crippen ? '--crippen' : ''}\
      ${params.threshold ? '--threshold ' + params.threshold : ''}

      # count the number of outputs - for some strange reason the fgrep command fails is the file is empty
      if [ -s 'o3da_${inputs.name}' ]
      then
        COUNT=\$(fgrep -c '\$\$\$\$' 'o3da_${inputs.name}')
      else
        COUNT=0
      fi
    """
}