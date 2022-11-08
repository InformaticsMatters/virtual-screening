params.publish_dir = ''
params.publish_dir_mode = 'move'
params.removehs = true
params.minimize_cycles = 500
params.rms_threshold = 1.0
params.interval = 100

process gen_conformers {

    container 'informaticsmatters/vs-moldb:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    path inputs

    output:
    path 'confs-*.cxsmi'

    """
    python -m moldb.conformers -i $inputs -o confs-${inputs.name}.cxsmi\
      --rms-threshold $params.rms_threshold\
      --minimize-cycles $params.minimize_cycles\
      --interval $params.interval
    """
}