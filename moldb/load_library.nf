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
    nextflow run moldb/load_library.nf --inputs data/10000.smi --library_name test123 --chunk_size 1000 -with-trace
*/

nextflow.enable.dsl=2


params.inputs = 'inputs.smi'
params.chunk_size = 10000

inputs = file(params.inputs)

// includes
include { split_txt } from '../nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { standardize } from '../nf-processes/rdkit/standardize.nf'
include { load_standardized } from '../nf-processes/moldb/db_load.nf'

// workflow definitions
workflow load_library {

    take:
    inputs

    main:
    split_txt(inputs)
    standardize(split_txt.out.flatten())
    load_standardized(standardize.out.collect())
}

workflow {
    load_library(inputs)
}
