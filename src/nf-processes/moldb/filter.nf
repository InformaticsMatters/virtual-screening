params.publish_dir = ''
params.publish_dir_mode = 'copy'
params.count = 10000


process extract_need_enum {

    container 'informaticsmatters/vs-moldb:2.0.0'
    publishDir params.publish_dir ?: '.', mode: params.publish_dir_mode, enabled: params.publish_dir as boolean

    input:
    path specification
    val outputfile // e.g. 'need-enum.smi'

    output:
    path "${outputfile}"

    script:
    """
    python -m moldb.filter --specification $specification --output-need-enum '$outputfile' --count $params.count
    """
}

process extract_need_conf {

    container 'informaticsmatters/vs-moldb:2.0.0'
    publishDir params.publish_dir ?: '.', mode: params.publish_dir_mode, enabled: params.publish_dir as boolean

    input:
    path specification
    val outputfile // e.g. 'need-enum.smi'

    output:
    path "${outputfile}"

    script:
    """
    python -m moldb.filter --specification $specification --output-need-conf '$outputfile' --count $params.count
    """
}