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


def prepare_lists(infile, outfile_enum, outfile_le_confs, data_dir):

    total = 0
    existing_enum = 0
    existing_le_confs = 0
    count_enum = 0
    count_le_confs = 0
    duplicates = 0
    errors = 0
    
    dups = set()
    
    utils.log_dm_event('Processing file', infile)
    with open(infile) as inf:
        utils.expand_path(outfile_enum)
        with open(outfile_enum, 'w') as outenum:
            utils.expand_path(outfile_le_confs)
            with open(outfile_le_confs, 'w') as outleconfs:

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

                    tgt_enum_smi = os.path.join(path, digest + '.smi')
                    tgt_enum_sdf = os.path.join(path, digest + '.sdf.gz')
                    tgt_le_confs = os.path.join(path, digest + '_le_confs.sdf.gz')

                    if smi in dups:
                        duplicates += 1
                        continue
                    else:
                        dups.add(smi)

                    if os.path.exists(tgt_enum_smi) and os.path.exists(tgt_enum_sdf):
                        existing_enum += 1
                    else:
                        outenum.write(smi + '\t' + uid + '\t' + digest + '\n')
                        count_enum += 1

                    if os.path.exists(tgt_le_confs):
                        existing_le_confs += 1
                    else:
                        outleconfs.write(smi + '\t' + uid + '\t' + digest + '\n')
                        count_le_confs += 1

    return total, existing_enum, existing_le_confs, count_enum, count_le_confs, duplicates, errors
                    
                    
def main():

    # Example:
    #   python3 prepare_enum_conf_lists.py -i foo.smi --outfile-enum need-enum.smi --outfile-le-confs need-confs.smi

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Prepare enumeration and conformer lists')
    parser.add_argument('-i', '--input', required=True, help="File with inputs")
    parser.add_argument('--outfile-enum', default='need-enum.smi', help="Output file for molecules needing enumeration")
    parser.add_argument('--outfile-confs', default='need-confs.smi',
                        help="Output file for molecules needing low energy 3D conformer generation")
    parser.add_argument('-d', '--data-dir', default='molecules/sha256', help="Directory with sharded data")

    args = parser.parse_args()
    utils.log("prepare_enum_conf_lists.py: ", args)

    total, existing_enum, existing_confs, count_enum, count_confs, duplicates, errors = \
        prepare_lists(args.input, args.outfile_enum, args.outfile_confs, args.data_dir)

    tmpl = 'Processed {} records. {} duplicates, {} errors.\
 {} already enumerated, {} already have low energy conformers,\
 {} need enumeration, {} need low energy conformers generated'

    utils.log_dm_event(tmpl.format(
        total, duplicates, errors, existing_enum, existing_confs, count_enum, count_confs))
    
    
if __name__ == "__main__":
    main()
        
