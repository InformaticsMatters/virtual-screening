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
Workflow that generates 3D alignments of ligands against fragment structures.

Example:
nextflow run frag-merge-pharmacophore.nf --inputs data/candidates.sdf\
  --fragments data/fragments.sdf --group_by_field std_smi --chunk_size 10\
  --publish_dir test
*/

nextflow.enable.dsl=2

/* params - other ones that can be specified are:
     split_sdf: chunk_size
     pharmacophore: rmsd, torsion_weight, count
     open3dalign: crippen, remove_hydrogens
*/
params.scratch = false
params.inputs = 'conformers.sdf'
params.fragments = 'fragments.sdf'
params.output_filename = 'ph4-results.sdf'
params.top = 100
params.gen3d = false
params.thresh_ph4 = null
params.thresh_o3da = null
params.group_by_field = null
params.publish_dir = './'

// includes
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { pharmacophore } from './nf-processes/plants/pharmacophore.nf'
include { open3dalign } from './nf-processes/rdkit/open3dalign.nf'
include { sucos } from './nf-processes/xchem/sucos.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf'
include { sd_best_sorted_top as filter } from './nf-processes/rdock/filter.nf'

def curr_t() {
    def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
    dateFormat.setTimeZone(TimeZone.getTimeZone('UTC'))
    return dateFormat.format(new java.util.Date())
}

// workflows
workflow ph4_align {

    take:
    inputs
    fragments

    main:
    def wrkflw = 'ph4_align'
    log.info("${curr_t()} # PROGRESS -START- $wrkflw:split_sdf 1")

    split_sdf(inputs)
    pharmacophore(split_sdf.out.flatten(), fragments, params.thresh_ph4 ?: 0)
    open3dalign(pharmacophore.out[0], fragments, params.thresh_o3da ?: 0)
    sucos(open3dalign.out[0], fragments)
    concatenate_files(sucos.out[0].collect(), params.output_filename, 'sucos_*.sdf')
    def filtered = channel.empty()
    if (params.group_by_field) {
        filtered = filter(concatenate_files.out, 'o3da_score_rel', true, params.group_by_field,
            params.output_filename[0..-5] + '-best.sdf')
    }

    def split_count = new java.util.concurrent.atomic.AtomicInteger()
    split_sdf.out.flatten().subscribe { _part ->
        def now = curr_t()
        if (split_count.get() == 0) log.info("$now # PROGRESS -DONE- $wrkflw:split_sdf 1")
        log.info("$now # PROGRESS -START- $wrkflw:pharmacophore ${split_count.incrementAndGet()}")
    }

    def pharmacophore_count = new java.util.concurrent.atomic.AtomicInteger()
    def open3dalign_count = new java.util.concurrent.atomic.AtomicInteger()
    def sucos_count = new java.util.concurrent.atomic.AtomicInteger()
    pharmacophore.out[1].subscribe { count_file ->
        def cost = count_file.text.trim() as Integer
        def now = curr_t()
        def n = pharmacophore_count.incrementAndGet()
        log.info("$now # PROGRESS -DONE- $wrkflw:pharmacophore $n")
        log.info("$now # PROGRESS -START- $wrkflw:open3dalign $n")
        log.info("$now # INFO -COST- +$cost ${n + open3dalign_count.get() + sucos_count.get()}")
    }

    open3dalign.out[1].subscribe { count_file ->
        def cost = count_file.text.trim() as Integer
        def now = curr_t()
        def n = open3dalign_count.incrementAndGet()
        log.info("$now # PROGRESS -DONE- $wrkflw:open3dalign $n")
        log.info("$now # INFO -COST- +$cost ${pharmacophore_count.get() + n + sucos_count.get()}")
    }

    sucos.out[1].subscribe { count_file ->
        def cost = count_file.text.trim() as Integer
        def now = curr_t()
        def n = sucos_count.incrementAndGet()
        log.info("$now # PROGRESS -DONE- $wrkflw:sucos $n")
        log.info("$now # INFO -COST- +$cost ${pharmacophore_count.get() + open3dalign_count.get() + n}")
    }

    sucos.out[0].collect().subscribe { _results ->
        log.info("${curr_t()} # PROGRESS -START- $wrkflw:concatenate_files 1")
    }

    concatenate_files.out.subscribe { _result ->
        def now = curr_t()
        log.info("$now # PROGRESS -DONE- $wrkflw:concatenate_files 1")
        if (params.group_by_field) {
            log.info("$now # PROGRESS -START- $wrkflw:filter 1")
        }
    }

    filtered.subscribe { _result ->
        log.info("${curr_t()} # PROGRESS -DONE- $wrkflw:filter 1")
    }

    emit:
    results = concatenate_files.out
    best = filtered
}

workflow {
    ph4_align(file(params.inputs), file(params.fragments))
}
