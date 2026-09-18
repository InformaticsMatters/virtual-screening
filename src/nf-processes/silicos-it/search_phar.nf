params.func_groups = '' // e.g. AROM,HDON,HACC,LIPO,CHARGE
params.threshold = 0
params.metric = 'TANIMOTO' // or TVERSKY_REF or TVERSKY_DB
params.merge = false
params.no_normal = false
params.no_hybrid = false
params.with_exclusion = false
params.epsilon = 0.5

process search_phar {

    // latest is the only tag they have. Not been updated for 7 years.
    container '3dechem/silicos-it:latest'

    input:
    path inputs // .sdf
    path query  // .mol or .phar

    output:
    path 'aligned_*.sdf'
    path 'aligned_*.tab'

    """
    align-it -r '$query' -d '$inputs'\
      ${query.name.endsWith('.phar') ? '--refType PHAR' : ''}\
      ${params.threshold ? '--cutOff ' + params.threshold + ' --rankBy \'' + params.metric + '\'' : ''}\
      --epsilon $params.epsilon\
      ${params.func_groups ? '--funcGroup ' + params.func_groups : ''}\
      ${params.merge ? '--merge' : ''}\
      ${params.no_normal ? --noNormal : ''}\
      ${params.no_hybrid ? --noHybrid : ''}\
      ${params.with_exclusion ? --withExclusion : ''}\
      --scores 'aligned_${inputs.name[0..-5]}.tab'\
      --out 'aligned_${inputs.name[0..-5]}.sdf'
    """
}


// process gen_phar {
//
//     container '3dechem/silicos-it:latest'
//
//     input:
//     path inputs
//
//     output:
//     tuple path inputs, path inputs.name[0..-5]}.phar
//
//     """
//     align-it -d '$inputs' -p '${inputs.name[0..-5]}.phar'
//     """
// }

// process search {
//
//     container '3dechem/silicos-it:latest'
//
//     input:
//     path ligand_phar
//     tuple path sdf, path phar
//
//     output:
//     path 'aligned_*.sdf'
//     path 'aligned_*.tab'
//
//     """
//     align-it -r '$ligand_phar' --refType PHAR -d '$sdf' -p '$phar'\
//           ${params.threshold ? '--cutOff ' + params.threshold + ' --rankBy \'' + params.metric + '\'': ''}\
//           --scores 'aligned_${sdf.name[0..-5]}.tab'\
//           --out 'aligned_${sdf.name[0..-5]}.sdf'
//     """
// }
