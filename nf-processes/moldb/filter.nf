params.publish_dir = ''
params.publish_dir_mode = 'copy'
params.output = 'need-enum.smi'
params.count = 10000


process extract_need_enum {

    container 'informaticsmatters/vs-moldb:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    file specification

    output:
    file params.output

    """
    python -m moldb.filter --specification $specification --output-need-enum $params.output --count $params.count
    """
}

process extract_need_conf {

    container 'informaticsmatters/vs-moldb:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    file specification

    output:
    file params.output

    """
    python -m moldb.filter --specification $specification --output-need-conf $params.output --count $params.count
    """
}