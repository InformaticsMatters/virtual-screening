params.publish_dir = ''
params.publish_dir_mode = 'copy'

/** Concatenate the files matching glob into outputfile.
*/
process concatenate_files {

    container 'informaticsmatters/vs-prep:2.0.0'
    publishDir params.publish_dir ?: '.', mode: params.publish_dir_mode, enabled: params.publish_dir as boolean

    input:
    path part
    val outputfile // e.g. 'results.sdf'
    val glob       // e.g. 'docked_*.sdf'

    output:
    path "${outputfile}"

    script:
    """
    DIR=\$(dirname "${outputfile}")
    mkdir -p \$DIR
    ls ${glob} | xargs cat >> ${outputfile}
    """
}
