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


import time, argparse
import utils, le_conformers
from dm_job_utilities.dm_log import DmLog

from rdkit import Chem


def execute(mol, outfile, rms_threshold=0.35, minimize_cycles=500, remove_hydrogens=False, num_conformers=None):

    utils.log('Executing ...')

    molconfs = le_conformers.gen_conformers(mol, rms_threshold, minimize_cycles, remove_hydrogens, num_conformers=num_conformers)
    with Chem.SDWriter(outfile) as writer:
        for idx in range(molconfs.GetNumConformers()):
            molconfs.SetDoubleProp('Energy', molconfs.GetConformer(idx).GetDoubleProp('Energy'))
            molconfs.SetDoubleProp('Energy_Delta', molconfs.GetConformer(idx).GetDoubleProp('Energy_Delta'))
            writer.write(molconfs, confId=idx)

    count = molconfs.GetNumConformers()
    utils.log('Generated', count, 'conformers')

    return count


### start main execution #########################################

def main():

    # Example:
    #  ./le_conformers_for_mol.py -i bar.mol -o conformers.sdf

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate conformers for molecule')
    parser.add_argument('-i', '--input', help="Input file as molfile")
    parser.add_argument('-s', '--smiles', help="Input file as SMILES string")
    parser.add_argument('-o', '--output', required=True, help="Output file as SD file")
    parser.add_argument('-n', '--num-conformers', type=int,
                        help="Number of conformers to generate. If not specified the Inhibox rules are used")
    parser.add_argument('-m', '--minimize-cycles', type=int, default=500, help="Number of MMFF minimisation cycles")
    parser.add_argument('-t', '--rms-threshold', type=float, default=1.0, help="RMS threshold for excluding conformers")
    parser.add_argument('--remove-hydrogens', action='store_true', help='Remove hydrogens from the output')

    args = parser.parse_args()
    DmLog.emit_event("le_conformers_for_mol: ", args)

    if args.input and args.smiles:
        raise ValueError("Can't specify both --input and --smiles options")

    if args.input:
        mol = Chem.MolFromMolFile(args.input)
    elif args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
    else:
        raise ValueError('Must specify input using one of the --input or --smiles options')

    if not mol:
        raise ValueError('Bad input molecule')

    start = time.time()
    count = execute(mol, args.output, minimize_cycles=args.minimize_cycles, num_conformers=args.num_conformers,
                    remove_hydrogens=args.remove_hydrogens, rms_threshold=args.rms_threshold)
    end = time.time()

    DmLog.emit_event('Generated:', count, 'conformers', 'Time (s):', end - start)
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()