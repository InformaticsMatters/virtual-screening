params.inputs = 'need-conf.smi'
params.data_dir = 'combined'
params.chunk_size = 1000
params.digits = 6
params.interval = 10000

inputsfile = file(params.inputs)
outputsdir = file(params.data_dir)

process splitter {

    container 'informaticsmatters/virt-screening-obabel:1.0.4'

    input:
    file inputs from inputsfile

    output:
    file 'x*.smi' into chunks

    """
    split -l $params.chunk_size -d -a $params.digits --additional-suffix .smi $inputs
    """
}

process gen_conformers {

    container 'informaticsmatters/virt-screening-obabel:1.0.4'

    input:
    file chunks from chunks.flatten()
    file data from outputsdir


    """
    python3 /data/work/jku/5V6S/gen_conformer.py -i $chunks --data-dir $data --interval $params.interval
    """
}

