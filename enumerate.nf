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


params.inputs = 'need-enum.smi'
params.data_dir = 'molecules/sha256'
params.chunk_size = 10000

// files
inputs_smi = file(params.inputs) // smiles with molecules to enumerate
data_dir = file(params.data_dir) // sharded data dir

// includes
include { split_txt } from './nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { enumerate } from './nf-processes/rdkit/enumerate.nf'

// workflow definitions
workflow enumerate_sharded {

    take:
    inputs_smi
    data_dir

    main:
    split_txt(inputs_smi)
    enumerate(split_txt.out.flatten(), data_dir)
}

workflow {
    enumerate_sharded(inputs_smi, data_dir)
}
