params.scratch = false
// docking params
params.num_dockings = 25
params.mode = 'dock' // or minimise, score etc. Look in /rDock_2013.1_src/data/scripts for the options


process rdock_docking {

    container 'informaticsmatters/vs-rdock:latest'
    errorStrategy 'ignore'
    maxRetries 3
    scratch params.scratch

    input:
    path part // e.g. mols_part_0019.sdf
    path protein
    path 'docking.prm'
    path 'docking.as'

    output:
    path 'rdock_part_*.sdf' optional true // e.g. rdock_part_0019.sdf
    path 'failed_part_*.sdf' optional true

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
      rbdock -r docking.prm -p '${params.mode}.prm' -n $params.num_dockings -i "\$f" -o "docked_\${f::-3}" > "rdock_out_\${f::-3}.log"
      if [ \$? != 0 ]; then
        cat "\$f" >> '${part.name.replace('ligands', 'failed')}'
      fi
    done

    # combine the results
    cat docked_mol*.sd > ${part.name.replace('mols_part_', 'rdock_part_')}
    """
}