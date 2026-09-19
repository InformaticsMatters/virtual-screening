params.publish_dir = ''
params.publish_dir_mode = 'copy'
params.interval = 10000


process calc_molprops {

    container 'informaticsmatters/vs-moldb:3.0.0'
    publishDir params.publish_dir ?: '.', mode: params.publish_dir_mode, enabled: params.publish_dir as boolean

    input:
    path inputs

    output:
    path 'calc_*'

    script:
    """
    python -m moldb.calc_molprops --input $inputs --output calc_$inputs --interval $params.interval
    """
}