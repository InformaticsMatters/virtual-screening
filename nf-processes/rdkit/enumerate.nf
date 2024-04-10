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

params.fragment = 'hac'
params.enumerate_charges = true
params.enumerate_chirals = true
params.enumerate_tautomers = true
params.combinatorial = false
params.interval = 100
params.try_embedding = true
params.add_hydrogens = false
params.max_tautomers = null
params.id_column = null
params.delimiter = null
params.min_hac = null
params.max_hac = null
params.min_ph = null // default is 5
params.max_ph = null // default is 9
params.min_charge = null
params.max_charge = null
params.num_charges = null

process enumerate {

    container 'informaticsmatters/vs-prep:stable'

    input:
    file inputs

    output:
    file "enumerated-*.sdf"
    env COUNT

    """
    python -m enumerate -i $inputs -o enumerated-${inputs.name}.sdf --interval $params.interval\
      ${params.enumerate_charges ? '--enumerate-charges' : ''}\
      ${params.enumerate_chirals ? '--enumerate-chirals' : ''}\
      ${params.enumerate_tautomers ? '--enumerate-tautomers' : ''}\
      ${params.combinatorial ? '--combinatorial' : ''}\
      ${params.max_tautomers ? '--max-tautomers ' + params.max_tautomers : ''}\
      ${params.fragment ? '--fragment \'' + params.fragment + '\'' : ''}\
      ${params.id_column ? '--id-column \'' + params.id_column + '\'' : ''}\
      ${params.delimiter ? '--delimiter \'' + params.delimiter + '\'' : ''}\
      ${params.min_hac ? '--min-hac ' + params.min_hac : ''}\
      ${params.max_hac ? '--max-hac ' + params.max_hac : ''}\
      ${params.min_ph ? '--min-ph ' + params.min_ph : ''}\
      ${params.max_ph ? '--max-ph ' + params.max_ph : ''}\
      ${params.min_charge ? '--min-charge ' + params.min_charge : ''}\
      ${params.max_charge ? '--max-charge ' + params.max_charge : ''}\
      ${params.num_charges ? '--num-charges ' + params.num_charges : ''}\
      ${params.try_embedding ? '--try-embedding' : ''}\
      ${params.add_hydrogens ? '--add-hydrogens' : ''}

    # count the number of outputs
    COUNT=\$(fgrep -c '\$\$\$\$' 'enumerated-${inputs.name}.sdf')
    """
}