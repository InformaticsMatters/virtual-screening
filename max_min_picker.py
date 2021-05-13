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


import argparse, time
import utils
from rdkit import Chem, SimDivFilters
from rdkit.Chem import rdMolDescriptors


def pick(input, output, count, threshold, interval=0):

    inputs = 0
    num_dups = 0
    
    data = []
    fingerprints = []
    
    duplicates = set()
    
    with open(input) as inf:
        with open(output, 'w') as outf:
            
            utils.log('Starting fingerprinting ...')
            t0 = time.time()
            for line in inf:
                inputs += 1
                
                if interval and inputs % interval == 0:
                    utils.log("... fingerprinted {} records".format(inputs))
                    
                tokens = line.strip().split('\t')
                
                smi = tokens[0]
                digest = tokens[1]
                if smi in duplicates:
                    num_dups += 1
                    continue
                else:
                    duplicates.add(smi)
                    
                mol = Chem.MolFromSmiles(smi)
                if not mol:
                    continue
                
                fingerprints.append(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol,2))
                data.append((smi, digest))
                
            t1 = time.time()
            utils.log('Fingerprinting took {} seconds'.format((t1-t0)))
            
            utils.log('Starting picking ...')
            t2 = time.time()
            mmp = SimDivFilters.MaxMinPicker()
            if not count:
                count = len(fingerprints)
            if threshold:
                picks, thresh = mmp.LazyBitVectorPickWithThreshold(fingerprints, len(fingerprints), count, 1.0 - threshold)
                utils.log('Final pick threshold was', 1.0 - thresh)
            else:
                picks = mmp.LazyBitVectorPick(fingerprints, len(fingerprints), count)
            
            t3 = time.time()
            utils.log('Picking took {} seconds'.format((t3-t2)))
            
            utils.log('Writing data ....')
            for pick in picks:
                d = data[pick]
                outf.write('{}\t{}\n'.format(d[0], d[1]))
            utils.log('Finished')
                
    return inputs, len(fingerprints), len(picks), num_dups



def main():

    # Example:
    #   python3 max_min_picker.py -i foo.smi -o bar.smi -c 1000

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Prepare enumeration and conformer lists')
    parser.add_argument('-i', '--input', required=True, help="File with inputs")
    parser.add_argument('-o', '--output', required=True, help="Output file")
    parser.add_argument('-c', '--count', type=int, help="Number to pick")
    parser.add_argument('-t', '--threshold', type=float, help="Similarity threshold")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    print("max_min_picker.py: ", args)
    
    if not args.count and not args.threshold:
        utils.log('Must specify count or threshold or both')
        exit(1)
    
    total, candidates, picked, dups = pick(args.input, args.output, args.count, args.threshold, interval=args.interval)
    print('Picked {} from {} molecules. {} duplicates'.format(picked, candidates, dups))
    
    
if __name__ == "__main__":
    main()
        
