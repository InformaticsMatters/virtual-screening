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


params.input = 'molecules.smi'
params.output = 'enumerated.sdf'

params.publish_dir = './'

/* Useful params of the modules
params.chunk_size = 1000 // chunk size for splitting
params.header = false    // smiles input has header line

params.try_embedding = true
params.add_hydrogens = false
params.max_tautomers = 25
params.id_column = null
params.min_hac = null
params.max_hac = null
params.min_ph = 5
params.max_ph = 9
params.min_charge = null
params.max_charge = null
params.num_charges = null
*/
// files
input = file(params.inputs) // smiles or sdf with molecules to enumerate

// includes
print("INPUT:" +  params.input)
if (input.name.endsWith('.sdf') || input.name.endsWith('.sdf.gz')) {
    include { split_sdf as splitter } from './nf-processes/file/split_sdf.nf'
} else {
    include { split_txt as splitter } from './nf-processes/file/split_txt.nf' addParams(suffix: '.smi')
}
include { enumerate } from './nf-processes/rdkit/enumerate.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.output,
    glob: 'enumerated-*.sdf')

dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))
def curr_t() { dateFormat.format(new java.util.Date()) }
int splits = 0

def now = curr_t()
def wrkflw = 'enumerate_mols'
log.info("$now # PROGRESS -START- $wrkflw:splitter 1")


// workflow definitions
workflow enumerate_mols {

    take:
    input

    main:

    splitter(input)
    enumerate(splitter.out.flatten())
    concatenate_files(enumerate.out[0].collect())

    split_count = 0
    splitter.out.flatten().subscribe {
        now = curr_t()
        if (split_count == 0) {
            log.info("$now # PROGRESS -DONE- $wrkflw:splitter 1")
        }
        split_count++
        log.info("$now # PROGRESS -START- $wrkflw:enumerate $split_count")
    }

    int cost = 0
    int enumerate_count = 0
    enumerate.out[1].subscribe {
        cost += new Integer(it)
        enumerate_count += 1
        now = curr_t()
        log.info("$now # INFO -COST- $cost $enumerate_count")
        log.info("$now # PROGRESS -DONE- $wrkflw:enumerate $enumerate_count")
        if (enumerate_count == split_count) {
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
    enumerate_mols(input)
}
