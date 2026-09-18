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
nextflow run rdock-docking.nf --ligands data/candidates.sdf --protein data/dhfr-receptor-ph7.mol2 \
     --prmfile data/docking.prm --asfile data/docking.as \
     --num_dockings 5 --chunk_size 10 --publish_dir ./test
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
    glob: 'docked_*.sdf')
include { concatenate_files as collect_failed } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output_basename + '_failed.sdf',
    glob: 'failed_*.sdf',
    optional: true)

/* COST events need to be formatted like this:
     2022-07-11T14:14:26+00:00 # INFO -COST- 10000 1
 PROGRESS events like this:
     2022-07-11T14:14:26+00:00 # PROGRESS -START- rdock_docking:rdock 10
     2022-07-11T14:14:26+00:00 # PROGRESS -DONE- rdock_docking:rdock 5
*/


dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))
def curr_t() { dateFormat.format(new java.util.Date()) }
int splits = 0

def now = curr_t()
def wrkflw = 'rdock_docking'
log.info("$now # PROGRESS -START- $wrkflw:split_sdf 1")


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

    split_sdf_count = 0
    split_sdf.out.flatten().subscribe {
        now = curr_t()
        if (split_sdf_count == 0) {
             log.info("$now # PROGRESS -DONE- $wrkflw:split_sdf 1")
        }
        split_sdf_count++
        log.info("$now # PROGRESS -START- $wrkflw:rdock $split_sdf_count")
    }

    int cost = 0
    int rdock_count = 0
    rdock.out[2].subscribe {
        cost += new Integer(it)
        rdock_count += 1
        now = curr_t()
        log.info("$now # INFO -COST- $cost $rdock_count")
        log.info("$now # PROGRESS -DONE- $wrkflw:rdock $rdock_count")
        log.info("$now # PROGRESS -START- $wrkflw:collect_results 1")
    }
    collect_results.out.subscribe {
        now = curr_t()
        log.info("$now # PROGRESS -DONE- $wrkflw:collect_results 1")
    }
    collect_failed.out.subscribe {
        now = curr_t()
        log.info("$now # PROGRESS -DONE- $wrkflw:collect_failed 1")
    }

    emit:
    collect_results.out
    collect_failed.out
}

workflow {
    rdock_docking(ligands_sdf, protein_mol2, prmfile, asfile)
}