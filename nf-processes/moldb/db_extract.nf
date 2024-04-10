params.count = 10000
params.outfile = 'outputs.smi'

process extract_molprops {

    container 'informaticsmatters/vs-moldb:stable'

    output:
    file params.outfile

    """
    python -m moldb.extract_need_molprops --outfile $params.outfile --count $params.count
    """
}