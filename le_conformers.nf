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

params.inputs = 'molecules.sdf'
params.output = 'conformers.sdf'

params.publish_dir = './'

/* Useful params of the modules
params.chunk_size = 1000 // chunk size for splitting
params.header = false    // smiles input has header line

params.removehs = false
params.minimize_cycles = 500
params.rms_threshold = 1.0
params.id_column = null
*/

inputs = file(params.inputs)

// includes
if (inputs.name.endsWith('.sdf') || inputs.name.endsWith('.sdf.gz')) {
    include { split_sdf as splitter } from './nf-processes/file/split_sdf.nf'
} else {
    include { split_txt as splitter } from './nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
}
include { gen_conformers } from './nf-processes/rdkit/gen_confs.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output,
    glob: 'confs-*.sdf')

dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))
def curr_t() { dateFormat.format(new java.util.Date()) }
int splits = 0

def now = curr_t()
def wrkflw = 'generate_conformers'
log.info("$now # PROGRESS -START- $wrkflw:splitter 1")


// workflow definitions
workflow generate_conformers {

    take:
    inputs

    main:
    splitter(inputs)
    gen_conformers(splitter.out.flatten())
    concatenate_files(gen_conformers.out[0].collect())

    split_count = 0
    splitter.out.flatten().subscribe {
        now = curr_t()
        if (split_count == 0) {
            log.info("$now # PROGRESS -DONE- $wrkflw:splitter 1")
        }
        split_count++
        log.info("$now # PROGRESS -START- $wrkflw:gen_conformers $split_count")
    }

    int cost = 0
    int count = 0
    gen_conformers.out[1].subscribe {
        count += 1
        cost += new Integer(it)
        now = curr_t()
        log.info("$now # INFO -COST- $cost $count")
        log.info("$now # PROGRESS -DONE- $wrkflw:gen_conformers $count")
        if (count == split_count) {
            log.info("$now # PROGRESS -START- $wrkflw:concatenate_files 1")
        }
    }
    concatenate_files.out.subscribe {
        now = curr_t()
        log.info("$now # PROGRESS -DONE- $wrkflw:concatenate_files 1")
    }

    emit:
    concatenate_files.out
}

workflow {
    generate_conformers(inputs)
}