/* Copyright 2024 Informatics Matters Ltd.

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

nextflow.enable.dsl=2

params.input = 'molecules.sdf'
params.output = 'conformers.sdf'

params.publish_dir = './'

/* Useful params of the modules
params.chunk_size = 1000 // chunk size for splitting
params.header = false    // smiles input has header line

params.num_conformers = null
params.removehs = false
params.minimize_cycles = 500
params.rms_threshold = 1.0
params.delimiter = null
params.id_column = null
*/

// includes
include { split_sdf } from './nf-processes/file/split_sdf.nf'
include { split_txt } from './nf-processes/file/split_txt.nf'
include { gen_conformers } from './nf-processes/rdkit/gen_confs.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf'

def curr_t() {
    def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
    dateFormat.setTimeZone(TimeZone.getTimeZone('UTC'))
    return dateFormat.format(new java.util.Date())
}


// workflow definitions
workflow generate_conformers {

    take:
    inputs

    main:
    def wrkflw = 'generate_conformers'
    log.info("${curr_t()} # PROGRESS -START- $wrkflw:splitter 1")

    def is_sdf = inputs.name.endsWith('.sdf') || inputs.name.endsWith('.sdf.gz')
    def parts = is_sdf ? split_sdf(inputs) : split_txt(inputs, '.smi')
    gen_conformers(parts.flatten())
    concatenate_files(gen_conformers.out[0].collect(), params.output, 'confs-*.sdf')

    def split_count = new java.util.concurrent.atomic.AtomicInteger()
    parts.flatten().subscribe { _part ->
        def now = curr_t()
        if (split_count.get() == 0) {
            log.info("$now # PROGRESS -DONE- $wrkflw:splitter 1")
        }
        log.info("$now # PROGRESS -START- $wrkflw:gen_conformers ${split_count.incrementAndGet()}")
    }

    def cost = new java.util.concurrent.atomic.AtomicInteger()
    def count = new java.util.concurrent.atomic.AtomicInteger()
    gen_conformers.out[1].subscribe { count_file ->
        def n = count.incrementAndGet()
        def total = cost.addAndGet(count_file.text.trim() as Integer)
        def now = curr_t()
        log.info("$now # INFO -COST- $total $n")
        log.info("$now # PROGRESS -DONE- $wrkflw:gen_conformers $n")
        if (n == split_count.get()) {
            log.info("$now # PROGRESS -START- $wrkflw:concatenate_files 1")
        }
    }
    concatenate_files.out.subscribe { _result ->
        log.info("${curr_t()} # PROGRESS -DONE- $wrkflw:concatenate_files 1")
    }

    emit:
    concatenate_files.out
}

workflow {
    generate_conformers(file(params.input))
}
