params.inputs = 'need-confs.smi'
params.data_dir = 'molecules/sha256'
params.removehs = true
params.minimize_cycles = 500
params.rms_threshold = 1.0
params.chunk_size = 100
params.digits = 6
params.interval = 10

inputsfile = file(params.inputs)
outputsdir = file(params.data_dir)

process splitter {

    container 'informaticsmatters/vs-prep:latest'

    input:
    file inputs from inputsfile

    output:
    file 'x*.smi' into chunks

    """
    split -l $params.chunk_size -d -a $params.digits --additional-suffix .smi $inputs
    """
}

process gen_conformers {

    container 'informaticsmatters/vs-prep:latest'

    input:
    file chunks from chunks.flatten()
    file data from outputsdir

    """
    /code/le_conformers.py -i $chunks --data-dir $data\
        ${params.removehs ? '--remove-hydrogens' : ''}\
        --rms-threshold $params.rms_threshold\
        --minimize-cycles $params.minimize_cycles\
        --interval $params.interval
    """
}
