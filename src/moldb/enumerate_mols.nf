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
params.specification = null

// outputs
params.file = 'need-enum.smi'

// filter options
// params.count = 10000

// split options
params.chunk_size = 1000

// includes
include { extract_need_enum } from '../nf-processes/moldb/filter.nf'
include { split_txt } from '../nf-processes/file/split_txt.nf'
include { enumerate } from '../nf-processes/moldb/enumerate.nf'
include { load_enum } from '../nf-processes/moldb/db_load.nf'

def curr_t() {
    def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
    dateFormat.setTimeZone(TimeZone.getTimeZone('UTC'))
    return dateFormat.format(new java.util.Date())
}

// workflow definitions
workflow enumerate_forms {

    take:
    specification

    main:
    def wrkflw = 'enumerate_forms'
    log.info("${curr_t()} # PROGRESS -START- $wrkflw:extract_need_enum 1")

    extract_need_enum(specification, params.file)
    split_txt(extract_need_enum.out, '.smi')
    enumerate(split_txt.out.flatten())
    load_enum(enumerate.out[0])

    extract_need_enum.out.subscribe { _extracted ->
        def now = curr_t()
        log.info("$now # PROGRESS -DONE- $wrkflw:extract_need_enum 1")
        log.info("$now # PROGRESS -START- $wrkflw:split_txt 1")
    }

    def split_count = new java.util.concurrent.atomic.AtomicInteger()
    split_txt.out.flatten().subscribe { _part ->
        def now = curr_t()
        if (split_count.get() == 0) log.info("$now # PROGRESS -DONE- $wrkflw:split_txt 1")
        log.info("$now # PROGRESS -START- $wrkflw:enumerate ${split_count.incrementAndGet()}")
    }

    def enumerate_count = new java.util.concurrent.atomic.AtomicInteger()
    enumerate.out[1].subscribe { _count_file ->
        def now = curr_t()
        def n = enumerate_count.incrementAndGet()
        log.info("$now # PROGRESS -DONE- $wrkflw:enumerate $n")
        log.info("$now # PROGRESS -START- $wrkflw:load_enum $n")
    }

    def load_enum_count = new java.util.concurrent.atomic.AtomicInteger()
    load_enum.out.subscribe { count_file ->
        def cost = count_file.text.trim() as Integer
        def n = load_enum_count.incrementAndGet()
        def now = curr_t()
        log.info("$now # PROGRESS -DONE- $wrkflw:load_enum $n")
        log.info("$now # INFO -COST- +$cost $n")
    }
}

workflow {
    enumerate_forms(file(params.specification))
}
