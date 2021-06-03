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


import argparse, os
import utils
from utils import get_path_from_digest


def prepare_3d_list(infile, outfile_enum, outfile_conf, data_dir):

    total = 0
    existing_enum = 0
    existing_conf = 0
    count_enum = 0
    count_conf = 0
    duplicates = 0
    errors = 0
    
    dups = set()
    
    with open(infile) as inf:
        with open(outfile_enum, 'w') as outenum:
            with open(outfile_conf, 'w') as outconf:
                for line in inf:
                    total += 1
                    tokens = line.strip().split('\t')
                    smi = tokens[0]
                    uid = tokens[1]
                    digest = tokens[2]
                    parts = [data_dir]
                    parts.extend(get_path_from_digest(digest))
                    path = os.path.join(*parts)
                    if not os.path.isdir(path):
                        utils.log('WARNING, path', path, 'not found')
                        errors += 1
                        continue
                    
                    tgt_enum = os.path.join(path, digest + '.smi')
                    tgt_conf = os.path.join(path, digest + '.sdf')
                    
                    if smi in dups:
                        duplicates += 1
                        continue
                    else:
                        dups.add(smi)
                        
                    if os.path.exists(tgt_enum):
                        existing_enum += 1
                    else:
                        outenum.write(smi + '\t' + uid + '\t' + digest + '\n')
                        count_enum += 1
                        
                    if os.path.exists(tgt_conf):
                        existing_conf += 1
                    else:
                        outconf.write(smi + '\t' + uid + '\t' + digest + '\n')
                        count_conf += 1
                
                
    return total, existing_enum, existing_conf, count_enum, count_conf, duplicates, errors
                    
                    
def main():

    # Example:
    #   python3 prepare_enum_conf_lists.py -i foo.smi --outfile-enum bar.smi --outfile-conf baz.smi -d combined

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Prepare enumeration and conformer lists')
    parser.add_argument('-i', '--input', required=True, help="File with inputs")
    parser.add_argument('--outfile-enum', required=True, help="Output file for molecules needing enumeration")
    parser.add_argument('--outfile-conf', required=True, help="Output file for molecules needing 3D conformer generation")
    parser.add_argument('-d', '--data-dir', required=True, help="Directory with data")

    args = parser.parse_args()
    utils.log("prepare_enum_conf_lists.py: ", args)
    
    total, existing_enum, existing_conf, count_enum, count_conf, duplicates, errors = prepare_3d_list(
        args.input, args.outfile_enum, args.outfile_conf, args.data_dir)
    utils.log('Processed {} records. {} already enumerated, {} already have conformer, {} needing enumeration, {} needing conformer generation, {} duplicates, {} errors'.format(
        total, existing_enum, existing_conf, count_enum, count_conf, duplicates, errors))
    
    
if __name__ == "__main__":
    main()
        
