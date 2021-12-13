#!/usr/bin/env nextflow

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