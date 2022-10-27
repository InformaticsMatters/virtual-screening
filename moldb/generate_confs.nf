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
   nextflow run moldb/generate_confs.nf -with-trace --min_hac 16 --max_hac 18 --min_rings 2 --min_aro_rings 1 --count 1000 --chunk_size 25
*/

nextflow.enable.dsl=2


params.file = 'need-confs.smi'
params.chunk_size = 1000

// filter options:
// params.count = 10000
// all the mol prop filters e.g. --min_hac 16

// includes
include { extract_need_conf } from '../nf-processes/moldb/filter.nf' addParams(output: params.file)
include { split_txt } from '../nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { gen_conformers } from '../nf-processes/moldb/gen_conformers.nf'
include { load_conf } from '../nf-processes/moldb/db_load.nf'

// workflow definitions
workflow gen_confs {

    main:
    extract_need_conf()
    split_txt(extract_need_conf.out)
    gen_conformers(split_txt.out.flatten())
    load_conf(gen_conformers.out)
}

workflow {
    gen_confs()
}
