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

params.inputs = 'need-confs.smi'
params.data_dir = 'molecules/sha256'

inputs_smi = file(params.inputs)
data_dir = file(params.data_dir)

// includes
include { split_txt } from './nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { gen_conformers } from './nf-processes/rdkit/gen_conformers.nf'

// workflow definitions
workflow generate_conformers {

    take:
    inputs_smi
    data_dir

    main:
    split_txt(inputs_smi)
    gen_conformers(split_txt.out.flatten(), data_dir)
}

workflow {
    generate_conformers(inputs_smi, data_dir)
}