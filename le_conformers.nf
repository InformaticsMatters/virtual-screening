/* Copyright 2022 Informatics Matters Ltd.

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

nextflow.enable.dsl=2

params.input = 'molecules.sdf'
params.output = 'conformers.sdf'

input = file(params.input)

// includes
include { split_txt } from './nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { gen_conformers } from './nf-processes/rdkit/gen_conformers.nf'

// workflow definitions
workflow generate_conformers {

    take:
    input

    main:
    split_txt(input)
    gen_conformers(split_txt.out.flatten())
}

workflow {
    generate_conformers(input)
}