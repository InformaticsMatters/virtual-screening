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
    nextflow run moldb/calc_molprops.nf --chunk_size 1000 -with-trace
*/

nextflow.enable.dsl=2


params.chunk_size = 10000
// other props used by modules
// params.count = 10000 // number of molecules to extract and calculate

// includes
include { extract_molprops } from '../nf-processes/moldb/db_extract.nf'
include { split_txt } from '../nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { calc_molprops } from '../nf-processes/moldb/calc_molprops.nf'
include { load_molprops } from '../nf-processes/moldb/db_load.nf'

def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))

def now = dateFormat.format(new java.util.Date())
def wrkflw = 'molprops'
log.info("$now # PROGRESS -START- $wrkflw:extract_molprops 1")

// workflow definitions
workflow molprops {

    main:
    extract_molprops()
    split_txt(extract_molprops.out)
    calc_molprops(split_txt.out.flatten())
    load_molprops(calc_molprops.out)

    extract_molprops.out.subscribe {
        now = dateFormat.format(new java.util.Date())
        log.info("$now # PROGRESS -DONE- $wrkflw:extract_molprops 1")
        log.info("$now # PROGRESS -START- $wrkflw:split_txt 1")
    }

    split_count = 0
    split_txt.out.flatten().subscribe {
        now = dateFormat.format(new java.util.Date())
        if (split_count == 0) log.info("$now # PROGRESS -DONE- $wrkflw:split_txt 1")
        split_count++
        log.info("$now # PROGRESS -START- $wrkflw:calc_molprops $split_count")
    }

    calc_count = 0
    calc_molprops.out.subscribe {
        now = dateFormat.format(new java.util.Date())
        calc_count++
        log.info("$now # PROGRESS -DONE- $wrkflw:calc_molprops $calc_count")
    }

    load_count = 0
    load_molprops.out.subscribe {
        cost = new Integer(it)
        load_count++
        now = dateFormat.format(new java.util.Date())
        log.info("$now # PROGRESS -DONE- $wrkflw:load_molprops $load_count")
        log.info("$now # INFO -COST- +$cost $load_count")
    }
}

workflow {
    molprops()
}
