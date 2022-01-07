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

import argparse, os, json
from hashlib import sha256
import utils

from rdkit import Chem, RDLogger

from standardize_molecule import standardize_to_iso_smiles

def process(inputs, data_dir, no_standardize=False, hits_file=None):
    num_hits = 0
    num_misses = 0

    if hits_file:
        hitsf = open(hits_file, "w")

    for smi in inputs:

        if no_standardize:
            std_smi = smi
        else:
            std_smi, mol = standardize_to_iso_smiles(smi)

        sha256_digest = sha256(bytes(std_smi, 'utf-8')).hexdigest()

        parts = [data_dir, 'sha256']
        parts.extend(utils.get_path_from_digest(sha256_digest))
        path = os.path.join(*parts)
        if not os.path.isdir(path):
            num_misses += 1
            continue
        molj = os.path.join(path, sha256_digest + '.json')
        if os.path.exists(molj):
            num_hits += 1
        else:
            num_misses += 1
            utils.log(smi)
            continue

        with open(molj) as f:
            data = json.load(f)
            uid = data['uuid']
            line = "{}\t{}\t{}".format(smi, uid, sha256_digest)
            utils.log(line)
            if hits_file:
                hitsf.write(line + '\n')

    if hits_file:
        hitsf.close()

    return num_hits, num_misses

def main():
    # Example usage:
    #  python find_molecules.py --smiles 'CCCC'

    RDLogger.logger().setLevel(RDLogger.ERROR)

    parser = argparse.ArgumentParser(description='Find molecules')
    parser.add_argument('-i', '--inputs', nargs='+', help='Input files')
    parser.add_argument('--hits-file', help='Output file for the hits')
    parser.add_argument('--data-dir', default='molecules', help="Data directory")
    parser.add_argument('--smiles', help="Smiles to search")
    parser.add_argument('--no-standardize', action='store_true', help='Molecules are already standardized')

    args = parser.parse_args()
    utils.log("find_molecules Args: ", args)


    if args.smiles:
        all_inputs = [ args.smiles ]

    all_inputs = []
    if args.inputs:
        for input in args.inputs:
            if input.endswith('.smi'):
                with open(input, 'r') as f:
                    lines = f.read().splitlines()
                    utils.log('Found', len(lines), 'inputs in', input)
                    all_inputs.extend(lines)
            elif input.endswith('.sdf'):
                supplr = Chem.ForwardSDMolSupplier(input)
                for mol in supplr:
                    all_inputs.append(Chem.MolToSmiles(mol))
            else:
                raise ValueError('Unsuppported file format: ' + input)

    # this does the processing
    num_hits, num_misses = process(all_inputs, args.data_dir,
                                   no_standardize=args.no_standardize, hits_file=args.hits_file)

    utils.log('Num hits:', num_hits, 'Num misses:', num_misses)


if __name__ == "__main__":
    main()
