params.publish_dir = ''
params.publish_dir_mode = 'copy'
params.glob = '*.sdf'
params.outputfile = 'results.sdf'
params.optional = false

process concatenate_files {

    container 'informaticsmatters/vs-prep:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    path part

    output:
    path params.outputfile optional params.optional

    """
    DIR=\$(dirname "${params.outputfile}")
    mkdir -p \$DIR
    ls ${params.glob} | xargs cat >> ${params.outputfile}
    """
}