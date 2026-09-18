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

// includes
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { rdock_docking as rdock } from './nf-processes/rdock/rdock_docking.nf'
include { concatenate_files as collect_results } from './nf-processes/file/concatenate_files.nf'
include { concatenate_files as collect_failed } from './nf-processes/file/concatenate_files.nf'

/* COST events need to be formatted like this:
     2022-07-11T14:14:26+00:00 # INFO -COST- 10000 1
 PROGRESS events like this:
     2022-07-11T14:14:26+00:00 # PROGRESS -START- rdock_docking:rdock 10
     2022-07-11T14:14:26+00:00 # PROGRESS -DONE- rdock_docking:rdock 5
*/

def curr_t() {
    def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
    dateFormat.setTimeZone(TimeZone.getTimeZone('UTC'))
    return dateFormat.format(new java.util.Date())
}


// workflows
workflow rdock_docking {

    take:
    ligands_sdf
    protein_mol2
    docking_prm
    docking_as

    main:
    def wrkflw = 'rdock_docking'
    log.info("${curr_t()} # PROGRESS -START- $wrkflw:split_sdf 1")

    split_sdf(ligands_sdf)
    rdock(split_sdf.out.flatten(), protein_mol2, docking_prm, docking_as)
    collect_results(rdock.out[0].collect(), params.output_basename + '.sdf', 'docked_*.sdf')
    collect_failed(rdock.out[1].collect(), params.output_basename + '_failed.sdf', 'failed_*.sdf')

    def split_sdf_count = new java.util.concurrent.atomic.AtomicInteger()
    split_sdf.out.flatten().subscribe { _part ->
        def now = curr_t()
        if (split_sdf_count.get() == 0) {
             log.info("$now # PROGRESS -DONE- $wrkflw:split_sdf 1")
        }
        log.info("$now # PROGRESS -START- $wrkflw:rdock ${split_sdf_count.incrementAndGet()}")
    }

    def cost = new java.util.concurrent.atomic.AtomicInteger()
    def rdock_count = new java.util.concurrent.atomic.AtomicInteger()
    rdock.out[2].subscribe { count_file ->
        def total = cost.addAndGet(count_file.text.trim() as Integer)
        def n = rdock_count.incrementAndGet()
        def now = curr_t()
        log.info("$now # INFO -COST- $total $n")
        log.info("$now # PROGRESS -DONE- $wrkflw:rdock $n")
        log.info("$now # PROGRESS -START- $wrkflw:collect_results 1")
    }
    collect_results.out.subscribe { _result ->
        log.info("${curr_t()} # PROGRESS -DONE- $wrkflw:collect_results 1")
    }
    collect_failed.out.subscribe { _result ->
        log.info("${curr_t()} # PROGRESS -DONE- $wrkflw:collect_failed 1")
    }

    emit:
    results = collect_results.out
    failed = collect_failed.out
}

workflow {
    rdock_docking(file(params.ligands), file(params.protein), file(params.prmfile), file(params.asfile))
}
