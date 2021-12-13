#!/usr/bin/env nextflow
/* Perform docking using smina, a fork of autodock vina.

Example:
nextflow run smina-docking.nf --ligands data/candidates.sdf --protein data/dhfr-receptor-ph7.pdbqt \
    --ligand data/dhfr-ligand.pdbqt --padding 2 --exhaustiveness 4 --scoring_function vinardo  \
    --publish_dir test --chunk_size 5
*/

nextflow.enable.dsl=2

params.chunk_size = 100
params.scratch = false

// docking params
params.ligands = 'ligands.sdf'
params.ligand = 'ligand.pdbqt'
params.protein = 'receptor.pdbqt'
params.outputfile = 'results_smina.sdf'


// files
ligands = file(params.ligands)
ligand = file(params.ligand)
protein = file(params.protein)


// includes
include { convert_format as format_protein } from './nf-processes/obabel/convert_format.nf' addParams(
    input_extensions: ['.pdb', '.mol2'],
    output_extension: '.pdbqt',
    outputfile: 'ready_receptor'
    )
include { convert_format as format_ligand } from './nf-processes/obabel/convert_format.nf' addParams(
    input_extensions: ['.pdb', '.mol2', '.mol'],
    output_extension: '.pdbqt',
    outputfile: 'ready_ligand'
    )
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { smina_docking as smina } from './nf-processes/smina/smina_docking.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf' addParams(
    glob: 'smina_*.sdf')


// workflows
workflow smina_docking {

    take:
    ligands
    ligand
    protein

    main:
    format_protein(protein)
    format_ligand(ligand)
    split_sdf(ligands)
    smina(split_sdf.out.flatten(), format_ligand.out, format_protein.out)
    concatenate_files(smina.out[0].collect())

    emit:
    concatenate_files.out
}

workflow {
    smina_docking(ligands, ligand, protein)
}

