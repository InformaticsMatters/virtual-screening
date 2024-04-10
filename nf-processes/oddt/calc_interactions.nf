params.scratch = false

params.key_hbond = false
params.key_hydrophobic = false
params.key_halogen = false
params.key_salt_bridge = false
params.key_pi_stacking = false
params.key_pi_cation = false
params.strict = true
params.exact_protein = false
params.exact_ligand = false

/* Generate a list of interactions using ODDT.
The inputs are expected to be a SD file (*.sdf) with the poses to analyse and a PDB file (*.pdb) with the protein
to which the poses correspond.
The output is a SD file corresponding to the inputs but with interaction data added as additional fields.
See the Python oddt_interactions.py module for full details.
*/
process calc_interactions {

    container 'informaticsmatters/vs-oddt:stable'
    errorStrategy 'retry'
    maxRetries 3
    scratch params.scratch

    input:
    path poses_sdf   // .sdf file
    path protein_pdb // .pdb file

    output:
    path 'oddt_*.sdf'

    """
    python /code/oddt_interactions.py -i '$poses_sdf' -p '$protein_pdb' -o 'oddt_${poses_sdf.name}'\
      ${params.key_hbond ? '--key-hbond \'' + params.key_hbond + '\'' : ''}\
      ${params.key_hydrophobic ? '--key-hydrophobic \'' + params.key_hydrophobic + '\'' : ''}\
      ${params.key_halogen ? '--key-halogen \'' + params.key_halogen + '\'' : ''}\
      ${params.key_salt_bridge ? '--key-salt-bridge \'' + params.key_salt_bridge + '\'' : ''}\
      ${params.key_pi_stacking ? '--key-pi-stacking \'' + params.key_pi_stacking + '\'' : ''}\
      ${params.key_pi_cation ? '--key-pi-cation \'' + params.key_pi_cation + '\'' : ''}\
      ${params.strict ? '--strict' : ''}\
      ${params.exact_protein ? '--exact-protein' : ''}\
      ${params.exact_ligand ? '--exact-ligand' : ''}
    """
}
