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
   nextflow run moldb/generate_confs.nf -with-trace --specification specification.txt --count 1000 --chunk_size 50
*/

nextflow.enable.dsl=2

// inputs
params.specification

params.file = 'need-confs.smi'
params.chunk_size = 1000

// filter options:
// params.count = 10000
// all the mol prop filters e.g. --min_hac 16

specification = file(params.specification)

// includes
include { extract_need_conf } from '../nf-processes/moldb/filter.nf' addParams(output: params.file)
include { split_txt } from '../nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
include { gen_conformers } from '../nf-processes/moldb/gen_conformers.nf'
include { load_conf } from '../nf-processes/moldb/db_load.nf'

def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))

def now = dateFormat.format(new java.util.Date())
def wrkflw = 'gen_confs'
log.info("$now # PROGRESS -START- $wrkflw:extract_need_conf 1")

// workflow definitions
workflow gen_confs {

    take:
    specification

    main:
    extract_need_conf(specification)
    split_txt(extract_need_conf.out)
    gen_conformers(split_txt.out.flatten())
    load_conf(gen_conformers.out)

    extract_need_conf.out.subscribe {
        now = dateFormat.format(new java.util.Date())
        log.info("$now # PROGRESS -DONE- $wrkflw:extract_need_conf 1")
        log.info("$now # PROGRESS -START- $wrkflw:split_txt 1")
    }

    split_count = 0
    split_txt.out.flatten().subscribe {
        now = dateFormat.format(new java.util.Date())
        if (split_count == 0) log.info("$now # PROGRESS -DONE- $wrkflw:split_txt 1")
        split_count++
        log.info("$now # PROGRESS -START- $wrkflw:gen_conformers $split_count")
    }

    conf_count = 0
    gen_conformers.out.subscribe {
        now = dateFormat.format(new java.util.Date())
        conf_count++
        log.info("$now # PROGRESS -DONE- $wrkflw:gen_conformers $conf_count")
    }

    load_count = 0
    load_conf.out.subscribe {
        cost = new Integer(it)
        load_count++
        now = dateFormat.format(new java.util.Date())
        log.info("$now # PROGRESS -DONE- $wrkflw:load_conf $load_count")
        log.info("$now # INFO -COST- +$cost $load_count")
    }
}

workflow {
    gen_confs(specification)
}
