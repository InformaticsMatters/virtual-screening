#!/usr/bin/env nextflow

params.chunk = 100
params.scratch = false

// docking params
params.ligands = 'ligands.sdf'
params.protein = 'receptor.mol2'
params.prmfile = 'docking.prm'
params.asfile = 'docking.as'
params.num_dockings = 25
params.publishDir = './'


// files
ligands = file(params.ligands)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)


process sdsplit {

    container 'informaticsmatters/vs-rdock:latest'

    input:
    file ligands

    output:
    file 'ligands_part*.sd' into ligand_parts

    """
    sdsplit -${params.chunk} -oligands_part_ $ligands

    for f in ligands_part_*.sd; do
      n=\${f:13:-3}
      if [ \${#n} == 1 ]; then
        mv \$f ligands_part_000\${n}.sd
      elif [ \${#n} == 2 ]; then
        mv \$f ligands_part_00\${n}.sd
      elif [ \${#n} == 3 ]; then
        mv \$f ligands_part_0\${n}.sd
      fi
    done
    """
}

process rdock {

    container 'informaticsmatters/vs-rdock:latest'
    errorStrategy 'ignore'
    maxRetries 3
    scratch params.scratch

    input:
    file part from ligand_parts.flatten()
    file protein
    file 'docking.prm' from prmfile
    file 'docking.as' from asfile

    output:
    file 'rdock_part_*.sd' optional true into docked_parts
    file 'failed_part_*.sd' optional true into failed_parts

    """
    # split into single molecules
    sdsplit -1 -omol $part

    # rename so that sorting is correct
    for f in mol*.sd; do
      n=\${f:3:-3}
      if [ \${#n} == 1 ]; then
        mv \$f mol00\${n}.sd
      elif [ \${#n} == 2 ]; then
        mv \$f mol0\${n}.sd
      fi
    done

    # do the docking
    for f in mol*.sd; do
      echo "Docking \$f"
      rbdock -r docking.prm -p dock.prm -n $params.num_dockings -i \$f -o docked_\${f::-3} > rdock_out_\${f::-3}.log
      if [ \$? != 0 ]; then
        cat \$f >> ${part.name.replace('ligands', 'failed')}
      fi
    done

    # combine the results
    cat docked_mol*.sd > ${part.name.replace('ligands', 'rdock')}
    """
}


process collect_and_report {

    container 'informaticsmatters/vs-rdock:latest'
    publishDir params.publishDir, mode: 'move'

    input:
    file part from docked_parts.collect()

    output:
    file 'results_rdock.sdf'
    
    """
    rm -f results_rdock.sdf
    ls rdock_*.sd | xargs cat >> results_rdock.sdf
    """
}


process collect_failures {

    container 'informaticsmatters/vs-rdock:latest'
    publishDir params.publishDir, mode: 'move'

    input:
    file failed from failed_parts.collect()

    output:
    file 'failed_rdock.sdf' optional true

    """
    rm -f failed_rdock.sdf
    ls failed_*.sd | xargs cat >> failed_rdock.sdf
    """
}

