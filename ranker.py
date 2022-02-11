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


import argparse, time

from rdkit import Chem

import utils


def rank_smi(filename, id_idx):
    count = 0
    data = {}
    with open(filename, 'r') as reader:
        for line in reader:
            count += 1
            tokens = line.split('\t')
            smiles = tokens[id_idx]
            data[smiles] = str(count)
    return data


def rank_sdf(filename, id_field):
    count = 0
    missing = 0
    data = {}
    with Chem.SDMolSupplier(filename) as supplr:
        for mol in supplr:
            count += 1
            if mol.HasProp(id_field):
                data[mol.GetProp(id_field)] = str(count)
            else:
                missing += 1
    return data, missing


def rank(inputs, output_tab, id_field):
    missing_total = 0
    all_data = None

    for input in inputs:
        if input.endswith('.sdf'):
            data, missing = rank_sdf(input, id_field)
            missing_total += missing
        elif input.endswith('.smi'):
            data = rank_smi(input, 0)
        if all_data is None: # first file
            all_data = {}
            for key in data:
                all_data[key] = [data[key]]
        else: # subsequent files
            for key in all_data:
                if key in data:
                    rank = data[key]
                else:
                    missing_total += 1
                    rank = None
                all_data[key].append(rank)

    for key in all_data:
        ranks = all_data[key]
        line = key + '\t' + '\t'.join(ranks)
        print(line)

    return len(all_data), missing_total


### start main execution #########################################

def main():

    # Example:
    #   ./ranker.py -i conformers.sdf --id-field _Name -o ranked.tab

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Rank scores')
    parser.add_argument('-i', '--inputs', nargs='+', required=True, help="Input molecules as SDF files")
    parser.add_argument('--id-field', required=True, help="Field name with IDs")
    parser.add_argument('-o', '--output', required=True, help="Output as tab separated file")

    args = parser.parse_args()
    utils.log_dm_event("ranker: ", args)

    start = time.time()
    count, missing = rank(args.inputs, args.output, args.id_field)
    end = time.time()

    utils.log_dm_event('Processed', count, 'records.', missing, 'missing. Time (s):', end - start)


if __name__ == "__main__":
    main()
