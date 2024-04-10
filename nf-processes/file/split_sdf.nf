/* Copyright 2024 Informatics Matters Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

params.chunk_size = 1000

process split_sdf {

    container 'informaticsmatters/vs-rdock:stable'

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
        mv \$f mols_part_0000\${n}.sdf
      elif [ \${#n} == 2 ]; then
        mv \$f mols_part_000\${n}.sdf
      elif [ \${#n} == 3 ]; then
        mv \$f mols_part_00\${n}.sdf
      elif [ \${#n} == 4 ]; then
        mv \$f mols_part_0\${n}.sdf
      else
        mv \$f mols_part_\${n}.sdf
      fi
    done
    """
}