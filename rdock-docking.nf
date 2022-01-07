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

/*
Run rDock docking on a set of candidate ligands provided in a SD file.

Example:
nextflow run rdock-docking.nf --ligands data/candidates.sdf --protein data/dhfr-receptor-ph7.mol2 --num_dockings 5\
    --chunk_size 10 --publish_dir ./test
*/

nextflow.enable.dsl=2

params.chunk_size = 100
params.scratch = false

// docking inputs
params.ligands = 'ligands.sdf'
params.protein = 'receptor.mol2'
params.prmfile = 'docking.prm'
params.asfile = 'docking.as'
params.output_basename = 'results_rdock'
params.publish_dir = './'

// files
ligands_sdf = file(params.ligands)
protein_mol2 = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)


// includes
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { rdock_docking as rdock } from './nf-processes/rdock/rdock_docking.nf'
include { concatenate_files as collect_results } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output_basename + '.sdf',
    glob: 'rdock_*.sdf')
include { concatenate_files as collect_failed } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output_basename + '_failed.sdf',
    glob: 'failed_*.sdf',
    optional: true)


// workflows
workflow rdock_docking {

    take:
    ligands_sdf
    protein_mol2
    docking_prm
    docking_as

    main:
    split_sdf(ligands_sdf)
    rdock(split_sdf.out.flatten(), protein_mol2, docking_prm, docking_as)
    collect_results(rdock.out[0].collect())
    collect_failed(rdock.out[1].collect())

    emit:
    collect_results.out
    collect_failed.out
}

workflow {
    rdock_docking(ligands_sdf, protein_mol2, prmfile, asfile)
}