params.publish_dir = ''
params.publish_dir_mode = 'copy'
params.interval = 10000


process calc_molprops {

    container 'informaticsmatters/vs-moldb:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    file inputs

    output:
    file 'calc_*'

    """
    python -m moldb.calc_molprops --input $inputs --output calc_$inputs --interval $params.interval
    """
}