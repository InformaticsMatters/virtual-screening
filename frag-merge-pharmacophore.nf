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
Workflow that generates 3D alignments of ligands against fragment structures.

Example:
nextflow run frag-merge-pharmacophore.nf --inputs data/candidates.sdf\
  --fragments data/Mpro-x0107_0A_hs.mol2 --group_by_field std_smi\
 --publish_dir test
*/

nextflow.enable.dsl=2

// params
// other ones that can be specified are: chunk_size, count, rmsd, torsion_weight
params.scratch = false
params.inputs = 'conformers.sdf'
params.fragments = 'fragments.sdf'
params.output_filename = 'ph4-results.sdf'
params.publish_dir = './'

// files
inputs = file(params.inputs)
fragments = file(params.fragments)

// includes
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { pharmacophore } from './nf-processes/plants/pharmacophore.nf'
include { open3dalign } from './nf-processes/rdkit/open3dalign.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output_filename,
    glob: 'o3da_*.sdf')
include { sd_best_sorted as filter } from './nf-processes/rdock/filter.nf' addParams(
    sort_field: 'o3da_score_rel',
    sort_descending: true,
    group_by_field: params.group_by_field,
    outputfile: params.output_filename[0..-5] + '-best.sdf'
)

// workflows
workflow ph4 {

    take:
    inputs
    fragments

    main:
    split_sdf(inputs)
    pharmacophore(split_sdf.out.flatten(), fragments)
    open3dalign(pharmacophore.out, fragments)
    concatenate_files(open3dalign.out.collect())
    if (params.group_by_field) {
        filter(concatenate_files.out)
    }

    emit:
    concatenate_files.out
    params.group_by_field ? filter.out : ''
}

workflow {
    ph4(inputs, fragments)
}