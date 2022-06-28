// process params
params.scratch = false
params.errorStrategy = 'retry'
params.maxRetries = 3

// docking params
params.num_dockings = 25
params.mode = 'dock' // or minimise, score etc. Look in /rDock_2013.1_src/data/scripts for the options


process rdock_docking {

    container 'informaticsmatters/vs-rdock:latest'
    //errorStrategy params.errorStrategy
    //maxRetries params.maxRetries
    errorStrategy { sleep(Math.pow(2, task.attempt) * 500 as long); return 'retry' }
    maxRetries 5
    scratch params.scratch

    input:
    path part    // e.g. mols_part_0019.sdf
    path protein // MOL2 format
    path 'docking.prm'
    path 'docking.as'

    output:
    path 'docked_*.sdf' optional true // e.g. docked_mols_part_0019.sdf
    path 'failed_*.sdf' optional true

    """
    set -e

    # split into single molecules
    sdsplit -1 -omol $part

    # rename so that sorting is correct
    for f in mol*.sd; do
      n=\${f:3:-3}
      if [ \${#n} == 1 ]; then
        mv \$f mol000\${n}.sd
      elif [ \${#n} == 2 ]; then
        mv \$f mol00\${n}.sd
      elif [ \${#n} == 3 ]; then
        mv \$f mol0\${n}.sd
      fi
    done

    # do the docking
    for f in mol*.sd; do
      echo "Docking \$f"
      rbdock -r docking.prm -p '${params.mode}.prm' -n $params.num_dockings -i "\$f" -o "rdock_\${f::-3}" > "rdock_out_\${f::-3}.log"
      if [ \$? != 0 ]; then
        cat "\$f" >> 'failed_${part.name}'
      fi
    done

    # combine the results
    cat rdock_*.sd > 'docked_${part.name}'
    """
}