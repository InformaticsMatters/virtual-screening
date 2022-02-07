#!/usr/bin/env python

# Copyright 2022 Informatics Matters Ltd.
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


import os, glob, argparse
import utils
from dm_job_utilities.dm_log import DmLog

def file_matches(filename, min_hac, max_hac):

    if min_hac is None and max_hac is None:
        return True
    sep_pos = filename.rindex(os.path.sep)
    num_part = int(filename[sep_pos+1:-4])
    if min_hac is not None and num_part < min_hac:
        return False
    if max_hac is not None and num_part > max_hac:
        return False
    return True
    

def list_files(input_dir, min_hac=None, max_hac=None):
    pattern = os.path.join(input_dir, '[0-9]*.smi')
    DmLog.emit_event('Inspecting files', pattern)
    files = glob.glob(pattern)
    matched = []
    for file in files:
        if file_matches(file, min_hac, max_hac):
            matched.append(file)
            
    return sorted(matched)


indexes = {}

def create_indexes(header):
    global indexes
    indexes = {}
    for i, token in enumerate(header):
        indexes[token] = i
        

def get_value(tokens, prop):
    idx = indexes[prop]
    return tokens[idx]


def check_value(tokens, prop, min_value, max_value):
    idx = indexes[prop]
    value = int(tokens[idx])
    if min_value is not None and min_value > value:
        return False
    if max_value is not None and max_value < value:
        return False
    return True
    
    

def filter(input_dirs, output_file, min_hac=None, max_hac=None,
    min_rotb=None, max_rotb=None,
    min_rings=None, max_rings=None,
    min_aro_rings=None, max_aro_rings=None,
    min_chiral_centres=None, max_chiral_centres=None,
    min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
    min_sp3=None, max_sp3=None):
    
    all_files = []
    for input_dir in input_dirs:
        files = list_files(input_dir, min_hac, max_hac)
        all_files.extend(files)
    DmLog.emit_event('Files to process:', all_files)
    count = 0
    total = 0
    num_dups = 0
    duplicates = set()

    utils.expand_path(output_file)
    with open(output_file, 'w') as out:
        for i, file in enumerate(all_files):
            DmLog.emit_event('Processing', file)
            with open(file) as f:
                header_line = f.readline()
                header_tokens = header_line.strip().split('\t')
                create_indexes(header_tokens)
                for line in f:
                    total += 1
                    tokens = line.strip().split('\t')
                    smi = get_value(tokens, 'smiles')
                    uid = get_value(tokens, 'uuid')
                    digest = get_value(tokens, 'sha256')
                    
                    if smi in duplicates:
                        num_dups += 1
                        continue
                    else:
                        duplicates.add(smi)
                    
                    if not check_value(tokens, 'rot_bonds', min_rotb, max_rotb):
                        continue
                        
                    if not check_value(tokens, 'ring_count', min_rings, max_rings):
                        continue
                        
                    if not check_value(tokens, 'aromatic_ring_count', min_aro_rings, max_aro_rings):
                        continue
                        
                    if not check_value(tokens, 'chiral_centres', min_chiral_centres, max_chiral_centres):
                        continue
                        
                    if not check_value(tokens, 'undefined_chiral_centres', min_undefined_chiral_centres, max_undefined_chiral_centres):
                        continue
                        
                    if not check_value(tokens, 'num_sp3', min_sp3, max_sp3):
                        continue
                    
                    count += 1
                    out.write('{}\t{}\t{}\n'.format(smi, uid, digest))
            
    return total, count, num_dups


def main():

    # Example:
    #   python3 -m filter -i chemspace_feb2021 -o foo.smi --min-hac 16 --max-hac 24 --min-rings 2 --min-aro-rings 1 --max-chiral-centres 2 --max-undefined-chiral-centres 0

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Filter')
    parser.add_argument('-i', '--inputs', nargs='+', required=True, help="Directories with inputs")
    parser.add_argument('-o', '--outfile', required=True, help="Output file")
    parser.add_argument('--min-hac', type=int, help="Min value for heavy atom count")
    parser.add_argument('--max-hac', type=int, help="Max value for heavy atom count")
    parser.add_argument('--min-rotb', type=int, help="Min value for rotatable bond count")
    parser.add_argument('--max-rotb', type=int, help="Max value for rotatable bond count")
    parser.add_argument('--min-rings', type=int, help="Min value for ring count")
    parser.add_argument('--max-rings', type=int, help="Max value for ring count")
    parser.add_argument('--min-aro-rings', type=int, help="Min value for aromatic ring count")
    parser.add_argument('--max-aro-rings', type=int, help="Max value for aromatic ring count")
    parser.add_argument('--min-chiral-centres', type=int, help="Min value for number of tetrahedral chiral centres")
    parser.add_argument('--max-chiral-centres', type=int, help="Max value for number of tetrahedral chiral centres")
    parser.add_argument('--min-undefined-chiral-centres', type=int,
                        help="Min value for number of undefined tetrahedral chiral centres")
    parser.add_argument('--max-undefined-chiral-centres', type=int,
                        help="Max value for number of undefined tetrahedral chiral centres")
    parser.add_argument('--min-sp3', type=int, help="Min value for SP3 count")
    parser.add_argument('--max-sp3', type=int, help="Max value for SP3 count")

    args = parser.parse_args()
    DmLog.emit_event("filter: ", args)
    
    total, count, num_dups = filter(args.inputs, args.outfile, 
        min_hac=args.min_hac, max_hac=args.max_hac,
        min_rotb=args.min_rotb, max_rotb=args.max_rotb,
        min_rings=args.min_rings, max_rings=args.max_rings,
        min_aro_rings=args.min_aro_rings, max_aro_rings=args.max_aro_rings,
        min_chiral_centres=args.min_chiral_centres, max_chiral_centres=args.max_chiral_centres,
        min_undefined_chiral_centres=args.min_undefined_chiral_centres, max_undefined_chiral_centres=args.max_undefined_chiral_centres,
        min_sp3=args.min_sp3, max_sp3=args.max_sp3)

    DmLog.emit_event('Matched {} out of {} records. {} duplicates'.format(count, total, num_dups))
    DmLog.emit_cost(count)
    

if __name__ == "__main__":
    main()
