params.publish_dir = ''
params.publish_dir_mode = 'move'
params.interval = 10000
params.num_charges = 2
params.try_embedding = true
params.add_hydrogens = false
params.max_tautomers = 25
params.output_format = 'cxsmi' // or sdf

process enumerate {

    container 'informaticsmatters/vs-prep:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    file inputs

    output:
    file "enumerated-*.${params.output_format}"

    """
    python -m moldb.enumerate -i $inputs -o enumerated-${inputs.name}.${params.output_format} --interval $params.interval\
      --enumerate-tautomers --enumerate-chirals --enumerate-charges\
      ${params.num_charges ? '--num-charges ' + params.num_charges : ''}\
      ${params.try_embedding ? '--try-embedding' : ''}\
      ${params.add_hydrogens ? '--add-hydrogens' : ''}\
      --max-tautomers $params.max_tautomers
    """
}