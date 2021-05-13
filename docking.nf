#!/usr/bin/env nextflow

params.scratch = false

// docking params
params.ligands = 'ligands/*.sdf'
params.protein = 'receptor.mol2'
params.prmfile = 'docking.prm'
params.asfile = 'docking.as'
params.num_dockings = 10


// interactions
params.iprotein = 'receptor.pdb'
params.key_hbond = null
params.key_hydrophobic = null
params.key_halogen = null
params.key_salt_bridge = null
params.key_pi_stacking = null
params.key_pi_cation = null


// files
ligands = file(params.ligands)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)

iprotein = file(params.iprotein)


process rdock {

    container 'informaticsmatters/rdock:2013.1'
    errorStrategy 'retry'
    maxRetries 3
    scratch params.scratch

    input:
    file part from ligands.flatten()
    file 'receptor.mol2' from protein
    file 'docking.prm' from prmfile
    file 'docking.as' from asfile

    output:
    file 'Docked_*.sd' optional true into docked_parts

    """
    rbdock -r $prmfile -p dock.prm -n $params.num_dockings -i $part -o ${part.name.replace('Prep', 'Docked')[0..-5]} > docked_out.log
    """
}


process interactions {

    container 'informaticsmatters/rdkit_pipelines:inters'
    errorStrategy 'retry'
    maxRetries 3
    scratch params.scratch

    input:
    file part from docked_parts
    file iprotein

    output:
    file 'INT_*.sdf' into interactions_parts

    """
    python -m pipelines.xchem.calc_interactions -i '$part' -if sdf -p $iprotein -o 'INT_${part.name[0..-4]}' -of sdf --no-gzip\
      ${params.key_hbond ? '--key-hbond ' + params.key_hbond : ''}\
      ${params.key_hydrophobic ? '--key-hydrophobic ' + params.key_hydrophobic : ''}\
      ${params.key_halogen ? '--key-halogen ' + params.key_halogen : ''}\
      ${params.key_salt_bridge ? '--key-salt-bridge ' + params.key_salt_bridge : ''}\
      ${params.key_pi_stacking ? '--key-pi-stacking ' + params.key_pi_stacking : ''}\
      ${params.key_pi_cation ? '--key-pi-cation ' + params.key_pi_cation : ''}\
      --nnscore /opt/python/NNScore_pdbbind2016.pickle\
      --rfscore /opt/python/RFScore_v1_pdbbind2016.pickle /opt/python/RFScore_v2_pdbbind2016.pickle /opt/python/RFScore_v3_pdbbind2016.pickle\
      --plecscore /opt/python/PLEClinear_p5_l1_pdbbind2016_s65536.pickle\
      --strict\
      --exact-ligand
    """
}

process filter_and_report {

    container 'informaticsmatters/rdock:2013.1'
    publishDir "./results", mode: 'copy'

    input:
    file part from interactions_parts.collect()

    output:
    file 'results.sdf.gz'
    file 'results.txt'

    """
    cat INT_*.sdf > results.sdf
    sdreport -t results.sdf > results.txt
    gzip results.sdf
    """
}

