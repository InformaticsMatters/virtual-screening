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

params.removehs = false
params.num_conformers = null
params.minimize_cycles = 500
params.rms_threshold = 1.0
params.delimiter = null
params.id_column = null
params.interval = 10

process gen_conformers {

    container 'informaticsmatters/vs-prep:stable'

    input:
    file inputs

    output:
    file "confs-*.sdf"
    env COUNT

    """
    /code/le_conformers.py -i $inputs -o confs-${inputs.name}.sdf\
        ${params.num_conformers ? '--num-conformers ' + params.num_conformers : ''}\
        ${params.removehs ? '--remove-hydrogens' : ''}\
        --rms-threshold $params.rms_threshold\
        --minimize-cycles $params.minimize_cycles\
        ${params.id_column ? '--id-column \'' + params.id_column + '\'' : ''}\
        ${params.delimiter ? '--delimiter \'' + params.delimiter + '\'' : ''}\
        --interval $params.interval

    # count the number of outputs
    COUNT=\$(fgrep -c '\$\$\$\$' 'confs-${inputs.name}.sdf')
    """
}