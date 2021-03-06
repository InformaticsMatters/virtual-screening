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

    """
    /code/open3dalign.py\
      --inputs $inputs\
      --query $query\
      --outfile o3da_${inputs.name}\
      --interval $params.interval\
      ${params.remove_hydrogens ? '--remove-hydrogens' : ''}\
      ${params.crippen ? '--crippen' : ''}\
      ${params.threshold ? '--threshold ' + threshold : ''}
    """
}