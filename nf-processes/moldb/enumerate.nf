params.enumerate_charges = true
params.enumerate_chirals = true
params.enumerate_tautomers = true
params.combinatorial = false
params.interval = 100
params.try_embedding = true
params.add_hydrogens = false
params.max_tautomers = null
params.min_ph = null // default is 5
params.max_ph = null // default is 9
params.min_charge = null
params.max_charge = null
params.num_charges = null


process enumerate {

    container 'informaticsmatters/vs-moldb:stable'

    input:
    file inputs

    output:
    file "enumerated-*.cxsmi"
    env COUNT

    """
    python -m moldb.enumerate -i $inputs -o enumerated-${inputs.name}.cxsmi --interval $params.interval\
      ${params.enumerate_charges ? '--enumerate-charges' : ''}\
      ${params.enumerate_chirals ? '--enumerate-chirals' : ''}\
      ${params.enumerate_tautomers ? '--enumerate-tautomers' : ''}\
      ${params.combinatorial ? '--combinatorial' : ''}\
      ${params.max_tautomers ? '--max-tautomers ' + params.max_tautomers : ''}\
      ${params.min_ph ? '--min-ph ' + params.min_ph : ''}\
      ${params.max_ph ? '--max-ph ' + params.max_ph : ''}\
      ${params.min_charge ? '--min-charge ' + params.min_charge : ''}\
      ${params.max_charge ? '--max-charge ' + params.max_charge : ''}\
      ${params.num_charges ? '--num-charges ' + params.num_charges : ''}\
      ${params.try_embedding ? '--try-embedding' : ''}\
      ${params.add_hydrogens ? '--add-hydrogens' : ''}

    # count the number of outputs
    COUNT=\$(wc -l < 'enumerated-${inputs.name}.cxsmi')
    """
}
