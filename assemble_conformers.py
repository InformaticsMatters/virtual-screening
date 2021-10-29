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


import argparse, os, time
import utils

from rdkit import Chem


def execute(input, output, data_dir, mode,
            exclude_base=False, exclude_tautomers=False, exclude_microstates=False, interval=0):

    if mode == 'single':
        ext = '.sdf'
    elif mode == 'low-energy':
        ext = '_le_confs.sdf'
    else:
        raise ValueError("Must specify a valid mode. 'single' or 'low-energy'")

    inputs = 0
    total = 0
    errors = 0
    duplicates = 0

    dups = set()

    with open(input) as inf:
        utils.expand_path(output)
        with open(output, 'w') as outf:
            for line in inf:
                inputs += 1

                if interval and inputs % interval == 0:
                    utils.log_dm_event("Processed {} records".format(inputs))

                tokens = line.strip().split('\t')
                smi = tokens[0]
                uid = tokens[1]
                digest = tokens[2]
                parts = [data_dir]
                parts.extend(utils.get_path_from_digest(digest))
                path = os.path.join(*parts)
                if not os.path.isdir(path):
                    utils.log_dm_event('WARNING, path', path, 'not found')
                    errors += 1
                    continue

                if smi in dups:
                    duplicates += 1
                    continue
                else:
                    dups.add(smi)

                total += 1

                confs = os.path.join(path, digest + ext)
                supplr = Chem.SDMolSupplier(confs)
                for mol in supplr:
                    if not mol:
                        errors += 1
                        continue
                    txt = supplr.GetItemText(0)
                    code = mol.GetProp('enum_code')
                    if exclude_base and code == 'B':
                        continue
                    if exclude_tautomers and (code == 'T' or code == 'X'):
                        continue
                    if exclude_microstates and (code == 'M' or code == 'X'):
                        continue
                    outf.write(txt)

    return inputs, total, errors, duplicates


def main():
    # Example:
    #   python3 assemble_conformers.py -i 20-25.smi -m single -o conformers.sdf

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Prepare enumeration and conformer lists')
    parser.add_argument('-i', '--input', required=True, help="File with inputs")
    parser.add_argument('-o', '--output', required=True, help="SDF file for outputs")
    parser.add_argument('-d', '--data-dir', default='molecules/sha256', help="Directory with sharded data")
    parser.add_argument('-m', '--mode', required=True, choices=['single', 'low-energy'],
                        help='Single conformer or low energy conformer mode [single, low-energy]')
    parser.add_argument('--exclude-base', action='store_true', help='Exclude base molecules')
    parser.add_argument('--exclude-tautomers', action='store_true', help='Exclude tautomers')
    parser.add_argument('--exclude-microstates', action='store_true', help='Exclude microstates')
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log_dm_event("prepare_enum_conf_lists.py: ", args)

    t0 = time.time()
    inputs, total, errors, duplicates = execute(args.input, args.output, args.data_dir, args.mode,
                                                exclude_base=args.exclude_base,
                                                exclude_tautomers=args.exclude_tautomers,
                                                exclude_microstates=args.exclude_microstates,
                                                interval=args.interval)
    t1 = time.time()

    utils.log_dm_event(
        'Processed {} inputs with {} unique mols in {}s. {} errors, {} duplicates'.format(inputs, total, (t1 - t0),
                                                                                          errors, duplicates))


if __name__ == "__main__":
    main()
