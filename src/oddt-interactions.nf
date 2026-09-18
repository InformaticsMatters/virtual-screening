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
Workflow that calculates interactions ligands make with a protein.
A workflow that takes a set of poses (.sdf file) and a protein (.pdb file) and calculates the interactions the
ligands make with the protein, and adds those interactions as extra properties to the output SD file.

Example:
nextflow run oddt-interactions.nf --ligands data/results_rdock.sdf.gz --protein data/dhfr-receptor.pdb --chunk_size 1000 --publish_dir test
*/

nextflow.enable.dsl=2

// params
params.scratch = false
params.ligands = 'ligands.sdf'
params.protein = 'receptor.pdb'
params.output_filename = 'oddt_interactions.sdf'
params.publish_dir = './'

// files
ligands = file(params.ligands)
protein = file(params.protein)

include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output_filename,
    glob: 'oddt_*.sdf')
include { calc_interactions } from './nf-processes/oddt/calc_interactions.nf'

workflow interactions {

    take:
    poses_sdf
    protein_pdb

    main:
    split_sdf(poses_sdf)
    calc_interactions(split_sdf.out.flatten(), protein_pdb)
    concatenate_files(calc_interactions.out.collect())

    emit:
    concatenate_files.out
}

workflow {
    interactions(ligands, protein)
}

