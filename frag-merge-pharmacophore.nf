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
     open3dalign: crippen, remove_hydrogens, threshold
     sd_best_sorted: group_by_field
*/
params.scratch = false
params.inputs = 'conformers.sdf'
params.fragments = 'fragments.sdf'
params.output_filename = 'ph4-results.sdf'
params.top = 100
params.gen3d = false
params.thresh_ph4 = null
params.thresh_o3da = null
params.publish_dir = './'

// files
inputs = file(params.inputs)
fragments = file(params.fragments)

// includes
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { pharmacophore } from './nf-processes/plants/pharmacophore.nf' addParams(
    threshold: params.thresh_ph4,
    gen3d: params.gen3d)
include { open3dalign } from './nf-processes/rdkit/open3dalign.nf' addParams(threshold: params.thresh_o3da)
include { concatenate_files } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output_filename,
    glob: 'o3da_*.sdf')
include { sd_best_sorted_top as filter } from './nf-processes/rdock/filter.nf' addParams(
    sort_field: 'o3da_score_rel',
    sort_descending: true,
    group_by_field: params.group_by_field,
    outputfile: params.output_filename[0..-5] + '-best.sdf'
)

def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))

def now = dateFormat.format(new java.util.Date())
def wrkflw = 'ph4_align'
log.info("$now # PROGRESS -START- $wrkflw:split_sdf 1")

// workflows
workflow ph4_align {

    take:
    inputs
    fragments

    main:
    split_sdf(inputs)
    pharmacophore(split_sdf.out.flatten(), fragments)
    open3dalign(pharmacophore.out[0], fragments)
    concatenate_files(open3dalign.out[0].collect())
    if (params.group_by_field) {
        filter(concatenate_files.out)
    }

    split_count = 0
    split_sdf.out.flatten().subscribe {
        now = dateFormat.format(new java.util.Date())
        if (split_count == 0) log.info("$now # PROGRESS -DONE- $wrkflw:split_sdf 1")
        split_count++
        log.info("$now # PROGRESS -START- $wrkflw:pharmacophore $split_count")
    }

    pharmacophore_count = 0
    open3dalign_count = 0
    pharmacophore.out[1].subscribe {
        cost = new Integer(it)
        now = dateFormat.format(new java.util.Date())
        pharmacophore_count++
        log.info("$now # PROGRESS -DONE- $wrkflw:pharmacophore $pharmacophore_count")
        log.info("$now # PROGRESS -START- $wrkflw:open3dalign $pharmacophore_count")
        log.info("$now # INFO -COST- +$cost ${pharmacophore_count + open3dalign_count}")
    }

    open3dalign.out[1].subscribe {
        cost = new Integer(it)
        now = dateFormat.format(new java.util.Date())
        open3dalign_count++
        log.info("$now # PROGRESS -DONE- $wrkflw:open3dalign $open3dalign_count")
        log.info("$now # INFO -COST- +$cost ${pharmacophore_count + open3dalign_count}")
    }

    open3dalign.out[0].collect().subscribe {
        now = dateFormat.format(new java.util.Date())
        log.info("$now # PROGRESS -START- $wrkflw:concatenate_files 1")
    }

    concatenate_files.out.subscribe {
        now = dateFormat.format(new java.util.Date())
        log.info("$now # PROGRESS -DONE- $wrkflw:concatenate_files 1")
        if (params.group_by_field) {
            log.info("$now # PROGRESS -START- $wrkflw:filter 1")
        }
    }

    if (params.group_by_field) {
        filter.out.subscribe {
            now = dateFormat.format(new java.util.Date())
            log.info("$now # PROGRESS -DONE- $wrkflw:filter 1")
        }
    }

    emit:
    concatenate_files.out
    params.group_by_field ? filter.out : ''
}

workflow {
    ph4_align(inputs, fragments)
}