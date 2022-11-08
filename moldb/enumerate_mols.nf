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

/* Example usage:
   nextflow run moldb/enumerate_mols.nf -with-trace --specification specification.txt --count 25000 --chunk_size 2500
*/

nextflow.enable.dsl=2

// inputs
params.specification

// outputs
params.file = 'need-enum.smi'

// filter options
// params.count = 10000

// split options
params.chunk_size = 1000

specification = file(params.specification)

// includes
include { extract_need_enum } from '../nf-processes/moldb/filter.nf' addParams(output: params.file)
include { split_txt } from '../nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { enumerate } from '../nf-processes/rdkit/enumerate.nf'
include { load_enum } from '../nf-processes/moldb/db_load.nf'

def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))

def now = dateFormat.format(new java.util.Date())
def wrkflw = 'enumerate_forms'
log.info("$now # PROGRESS -START- $wrkflw:extract_need_enum 1")

// workflow definitions
workflow enumerate_forms {

    take:
    specification

    main:
    extract_need_enum(specification)
    split_txt(extract_need_enum.out)
    enumerate(split_txt.out.flatten())
    load_enum(enumerate.out)

    extract_need_enum.out.subscribe {
        now = dateFormat.format(new java.util.Date())
        log.info("$now # PROGRESS -DONE- $wrkflw:extract_need_enum 1")
        log.info("$now # PROGRESS -START- $wrkflw:split_txt 1")
    }

    split_count = 0
    split_txt.out.flatten().subscribe {
        now = dateFormat.format(new java.util.Date())
        if (split_count == 0) log.info("$now # PROGRESS -DONE- $wrkflw:split_txt 1")
        split_count++
        log.info("$now # PROGRESS -START- $wrkflw:enumerate $split_count")
    }

    enumerate_count = 0
    enumerate.out.subscribe {
        now = dateFormat.format(new java.util.Date())
        enumerate_count++
        log.info("$now # PROGRESS -DONE- $wrkflw:enumerate $enumerate_count")
    }

    load_enum_count = 0
    load_enum.out.subscribe {
        cost = new Integer(it)
        load_enum_count++
        now = dateFormat.format(new java.util.Date())
        log.info("$now # PROGRESS -DONE- $wrkflw:load_enum $load_enum_count")
        log.info("$now # INFO -COST- +$cost $load_enum_count")
    }
}

workflow {
    enumerate_forms(specification)
}
