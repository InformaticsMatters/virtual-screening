#!/usr/bin/env nextflow

params.chunk = 25
params.scratch = false

// docking params
params.ligands = 'ligands.sdf'
params.protein = 'receptor.mol2'
params.prmfile = 'docking.prm'
params.asfile = 'docking.as'
params.num_dockings = 50


// files
ligands = file(params.ligands)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)


process sdsplit {

    container 'informaticsmatters/rdock:2013.1'

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

    container 'informaticsmatters/rdock:2013.1'
    errorStrategy 'retry'
    maxRetries 3
    scratch params.scratch

    input:
    file part from ligand_parts.flatten()
    file 'receptor.mol2' from protein
    file 'docking.prm' from prmfile
    file 'docking.as' from asfile

    output:
    file 'Docked_*.sd' optional true into docked_parts

    """
    rbdock -r $prmfile -p dock.prm -n $params.num_dockings -i $part -o ${part.name.replace('ligands', 'Docked')[0..-4]} > docked_out.log
    """
}


process collect_and_report {

    container 'informaticsmatters/rdock:2013.1'
    publishDir "./results", mode: 'move'

    input:
    file part from docked_parts.collect()

    output:
    file 'results.sdf.gz'
    file 'results.txt'

    """
    rm -f results.sdf
    ls Docked_*.sd | xargs cat >> results.sdf
    sdreport -t results.sdf > results.txt
    gzip results.sdf
    """
}

