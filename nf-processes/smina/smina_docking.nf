params.scratch = false

params.retries = 3
params.scoring_function = 'vina'
params.exhaustiveness = 8
params.padding = 4
params.cpu = 1

process smina_docking {

    container 'informaticsmatters/vs-smina:stable'
    errorStrategy 'retry'
    maxRetries params.retries
    scratch params.scratch

    input:
    file part
    file 'ligand.pdbqt'
    file 'receptor.pdbqt'

    output:
    file 'smina_part_*.sdf'

    """
    smina -r receptor.pdbqt -l $part --autobox_ligand ligand.pdbqt --autobox_add $params.padding\
      --exhaustiveness $params.exhaustiveness --scoring $params.scoring_function --cpu $params.cpu\
      -o ${part.name.replace('mols_', 'smina_')} > smina_out.log
    """
}