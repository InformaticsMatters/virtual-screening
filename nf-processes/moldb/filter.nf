params.publish_dir = ''
params.publish_dir_mode = 'move'
params.output = 'need-enum.smi'
params.count = 10000

params.min_hac = null
params.max_hac = null
params.min_rotb = null
params.max_rotb = null
params.min_rings = null
params.max_rings = null
params.min_aro_rings = null
params.max_aro_rings = null
params.min_chiral_centres = null
params.max_chiral_centres = null
params.min_undefined_chiral_centres = null
params.max_undefined_chiral_centres = null
params.min_sp3 = null
params.max_sp3 = null
params.min_logp = null
params.max_logp = null
params.min_tpsa = null
params.max_tpsa = null


process extract_need_enum {

    container 'informaticsmatters/vs-moldb:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    output:
    file params.output

    """
    python -m moldb.filter --output-need-enum $params.output --count $params.count\
    ${params.min_hac ? '--min-hac ' + params.min_hac : ''}\
    ${params.max_hac ? '--max-hac ' + params.max_hac : ''}\
    ${params.min_rotb ? '--min-rotb ' + params.min_rotb : ''}\
    ${params.max_rotb ? '--max-rotb ' + params.max_rotb : ''}\
    ${params.min_rings ? '--min-rings ' + params.min_rings : ''}\
    ${params.max_rings ? '--max-rings ' + params.max_rings : ''}\
    ${params.min_aro_rings ? '--min-aro-rings ' + params.min_aro_rings : ''}\
    ${params.max_aro_rings ? '--max-aro-rings ' + params.max_aro_rings : ''}\
    ${params.min_chiral_centres ? '--min-chiral-centres ' + params.min_chiral_centres : ''}\
    ${params.max_chiral_centres ? '--max-chiral-centres ' + params.max_chiral_centres : ''}\
    ${params.min_undefined_chiral_centres ? '--min-undefined-chiral-centres ' + params.min_undefined_chiral_centres : ''}\
    ${params.max_undefined_chiral_centres ? '--max-undefined-chiral-centres ' + params.max_undefined_chiral_centres : ''}\
    ${params.min_sp3 ? '--min-sp3 ' + params.min_sp3 : ''}\
    ${params.max_sp3 ? '--max-sp3 ' + params.max_sp3 : ''}\
    ${params.min_logp ? '--min-logp ' + params.min_logp : ''}\
    ${params.max_logp ? '--max-logp ' + params.max_logp : ''}\
    ${params.min_tpsa ? '--min-tpsa ' + params.min_tpsa : ''}\
    ${params.max_tpsa ? '--max-tpsa ' + params.max_tpsa : ''}
    """
}

process extract_need_conf {

    container 'informaticsmatters/vs-moldb:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    output:
    file params.output

    """
    python -m moldb.filter --output-need-conf $params.output --count $params.count\
    ${params.min_hac ? '--min-hac ' + params.min_hac : ''}\
    ${params.max_hac ? '--max-hac ' + params.max_hac : ''}\
    ${params.min_rotb ? '--min-rotb ' + params.min_rotb : ''}\
    ${params.max_rotb ? '--max-rotb ' + params.max_rotb : ''}\
    ${params.min_rings ? '--min-rings ' + params.min_rings : ''}\
    ${params.max_rings ? '--max-rings ' + params.max_rings : ''}\
    ${params.min_aro_rings ? '--min-aro-rings ' + params.min_aro_rings : ''}\
    ${params.max_aro_rings ? '--max-aro-rings ' + params.max_aro_rings : ''}\
    ${params.min_chiral_centres ? '--min-chiral-centres ' + params.min_chiral_centres : ''}\
    ${params.max_chiral_centres ? '--max-chiral-centres ' + params.max_chiral_centres : ''}\
    ${params.min_undefined_chiral_centres ? '--min-undefined-chiral-centres ' + params.min_undefined_chiral_centres : ''}\
    ${params.max_undefined_chiral_centres ? '--max-undefined-chiral-centres ' + params.max_undefined_chiral_centres : ''}\
    ${params.min_sp3 ? '--min-sp3 ' + params.min_sp3 : ''}\
    ${params.max_sp3 ? '--max-sp3 ' + params.max_sp3 : ''}\
    ${params.min_logp ? '--min-logp ' + params.min_logp : ''}\
    ${params.max_logp ? '--max-logp ' + params.max_logp : ''}\
    ${params.min_tpsa ? '--min-tpsa ' + params.min_tpsa : ''}\
    ${params.max_tpsa ? '--max-tpsa ' + params.max_tpsa : ''}
    """
}