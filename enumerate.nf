params.inputs = 'inputs.smi'
params.data_dir = 'combined'
params.chunk_size = 10000
params.digits = 6
params.interval = 1000

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

process enumerate {

    container 'informaticsmatters/vs-prep:latest'

    input:
    file chunks from chunks.flatten()
    file data from outputsdir


    """
    /code/enumerate.py -i $chunks --data-dir $data --interval $params.interval --enumerate-tautomers --enumerate-chirals --enumerate-charges
    """
}

