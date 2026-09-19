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

/**
Workflow that generates 3D alignments of ligands against a reference structure.

Example:
nextflow run open3dalign.nf --inputs data/candidates.sdf --query data/dhfr-lignad.mol --publish_dir test
*/

nextflow.enable.dsl=2

// params
// other ones that can be specified are: chunk_size, crippen, remove_hydrogens
params.scratch = false
params.inputs = 'conformers.sdf'
params.query = 'ligand.mol'
params.output_filename = 'open3dalign.sdf'
params.publish_dir = './'
params.group_by_field = null
params.threshold = 0

// includes
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { open3dalign } from './nf-processes/rdkit/open3dalign.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf'
include { sd_best_sorted as filter } from './nf-processes/rdock/filter.nf'

// workflows
workflow o3da {

    take:
    inputs
    query

    main:
    split_sdf(inputs)
    open3dalign(split_sdf.out.flatten(), query, params.threshold ?: 0)
    concatenate_files(open3dalign.out[0].collect(), params.output_filename, 'o3da_*.sdf')
    def filtered = channel.empty()
    if (params.group_by_field) {
        filtered = filter(concatenate_files.out, 'o3da_score_rel', true, params.group_by_field,
            params.output_filename[0..-5] + '-best.sdf')
    }

    emit:
    results = concatenate_files.out
    best = filtered
}

workflow {
    o3da(file(params.inputs), file(params.query))
}
