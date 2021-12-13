#!/usr/bin/env nextflow

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
