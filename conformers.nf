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

/* Example usage:
   nextflow run conformers.nf -with-trace -with-report --publish_dir outputs/conformers
*/

nextflow.enable.dsl=2


params.inputs = 'need-conf.smi'
params.chunk_size = 10000

// files
inputs_smi = file(params.inputs) // smiles with molecules to enumerate

// includes
include { split_txt } from './nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { enumerate } from './nf-processes/rdkit/enumerate.nf'

// workflow definitions
workflow enumerate_forms {

    take:
    inputs_smi

    main:
    split_txt(inputs_smi)
    enumerate(split_txt.out.flatten())

    emit:
    enumerate.out
}

workflow {
    enumerate_forms(inputs_smi)
}
