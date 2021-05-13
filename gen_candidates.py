#!/usr/bin/env python

# Copyright 2021 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import argparse, os, shutil, time
import utils


def execute(input, output, data_dir, interval=0):
    inputs = 0
    total = 0
    errors = 0
    duplicates = 0
    
    dups = set()
    
    with open(input) as inf:
        with open(output, 'w') as outf:
            for line in inf:
                inputs += 1
                
                if interval and inputs % interval == 0:
                    utils.log("Processed {} records".format(inputs))
                
                tokens = line.strip().split('\t')
                smi = tokens[0]
                digest = tokens[1]
                parts = [data_dir]
                parts.extend(utils.get_path_from_digest(digest))
                path = os.path.join(*parts)
                if not os.path.isdir(path):
                    utils.log('WARNING, path', path, 'not found')
                    errors += 1
                    continue
                    
                if smi in dups:
                    duplicates += 1
                    continue
                else:
                    dups.add(smi)
                
                total += 1
                    
                conf = os.path.join(path, digest + '.sdf')
                with open(conf) as sdf:
                    shutil.copyfileobj(sdf, outf)
                
    return inputs, total, errors, duplicates
    

def main():

    # Example:
    #   python3 gen_candidates.py -i 20-25.smi -0 candidates.sdf -d combined

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Prepare enumeration and conformer lists')
    parser.add_argument('-i', '--input', required=True, help="File with inputs")
    parser.add_argument('-o', '--output', required=True, help="SDF file for outputs")
    parser.add_argument('-d', '--data-dir', required=True, help="Directory with data")
    parser.add_argument("--interval", type=int, help="Reporting interval")


    args = parser.parse_args()
    utils.log("prepare_enum_conf_lists.py: ", args)
    
    t0 = time.time()
    inputs, total, errors, duplicates  = execute(args.input, args.output, args.data_dir, interval=args.interval)
    t1 = time.time()
    
    utils.log('Processed {} inputs with {} unique mols in {}s. {} errors, {} duplicates'.format(inputs, total, (t1 - t0), errors, duplicates))
    
    
if __name__ == "__main__":
    main()
