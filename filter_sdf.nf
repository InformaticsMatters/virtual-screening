nextflow.enable.dsl=2

params.inputs = 'inputs.sdf'

inputs = file(params.inputs)

// includes
include { sd_best_sorted as filter } from './nf-processes/rdock/filter.nf'

workflow filter_sdf {
    take:
    inputs_sdf

    main:
    filter(inputs_sdf)

}

workflow {
    filter_sdf(inputs)
}