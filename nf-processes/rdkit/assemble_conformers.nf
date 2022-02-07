
params.mode = 'single'
params.input = 'inputs.smi'
params.output = 'conformers.sdf'
params.interval = 10000

process assemble {

    container 'informaticsmatters/vs-prep:latest'

    input:
    path input // .smi
    path data_dir  // e.g. molecules/sha256

    output:
    path params.output

    """
    /code/assemble_conformers.py -m $params.mode\
      --input $input\
      --data-dir $data_dir\
      --output $params.output
    """
}