params.chunk_size = 5000

process split_sdf {

    container 'informaticsmatters/vs-rdock:latest'

    input:
    path molecules

    output:
    path 'mols_part*.sdf'

    """
    file=$molecules.name
    if [ \${file##*.} == 'gz' ]; then
        zcat $molecules | sdsplit -${params.chunk_size} -omols_part_
    else
        sdsplit -${params.chunk_size} -omols_part_ $molecules
    fi

    for f in mols_part_*.sd; do
      n=\${f:10:-3}
      if [ \${#n} == 1 ]; then
        mv \$f mols_part_000\${n}.sdf
      elif [ \${#n} == 2 ]; then
        mv \$f mols_part_00\${n}.sdf
      elif [ \${#n} == 3 ]; then
        mv \$f mols_part_0\${n}.sdf
      fi
    done
    """
}