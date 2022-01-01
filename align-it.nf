/** Workflow for pharmacophore based alignment using Silicos-it's align-it tool.
See http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/align-it/1.0.4/align-it.html for more details
on align-it.
This workflow is NOT performance optimised as the pharmacophores are generated as part of the search.
*/

nextflow.enable.dsl=2

params.query = 'ligand.mol' // or can be a .phar file
params.input = 'conformers.sdf'
params.output_scores = 'results.tab'
params.output_sdf = 'results.sdf'
params.chunk_size = 5000
params.publish_dir = './'

query = file(params.query)
inputs = file(params.input)


// includes
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { search_phar } from './nf-processes/silicos-it/search_phar.nf'
include { concatenate_files as concat_sdf } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output_sdf,
    glob: 'aligned_*.sdf')
include { concatenate_files as concat_scores } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output_scores,
    glob: 'aligned_*.tab')


// workflows
workflow align_it {

    take:
    query_mol_or_phar
    inputs_sdf

    main:
    split_sdf(inputs_sdf)
    search_phar(split_sdf.out.flatten(), query_mol_or_phar)
    concat_sdf(search_phar.out[0].collect())
    concat_scores(search_phar.out[1].collect())

    emit:
    concat_sdf.out
    concat_scores.out
}

workflow {
    align_it(query, inputs)
}
