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