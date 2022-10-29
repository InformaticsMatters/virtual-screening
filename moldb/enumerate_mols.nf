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
   nextflow run moldb/enumerate_mols.nf -with-trace --min_hac 16 --max_hac 20 --min_rings 2 --min_aro_rings 1
*/

nextflow.enable.dsl=2

// filter options
params.file = 'need-enum.smi'
// params.count = 10000
// all the mol prop filters e.g. --min_hac 16

// split options
params.chunk_size = 1000


// includes
include { extract_need_enum } from '../nf-processes/moldb/filter.nf' addParams(output: params.file)
include { split_txt } from '../nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { enumerate } from '../nf-processes/rdkit/enumerate.nf'
include { load_enum } from '../nf-processes/moldb/db_load.nf'

// workflow definitions
workflow enumerate_forms {

    main:
    extract_need_enum()
    split_txt(extract_need_enum.out)
    enumerate(split_txt.out.flatten())
    load_enum(enumerate.out)
}

workflow {
    enumerate_forms()
}
