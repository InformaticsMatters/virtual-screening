params.removehs = true
params.minimize_cycles = 500
params.rms_threshold = 1.0
params.interval = 10

process gen_conformers {

    container 'informaticsmatters/vs-prep:latest'

    input:
    file inputs
    file data_dir

    """
    /code/le_conformers.py -i $inputs --data-dir $data_dir\
        ${params.removehs ? '--remove-hydrogens' : ''}\
        --rms-threshold $params.rms_threshold\
        --minimize-cycles $params.minimize_cycles\
        --interval $params.interval
    """
}