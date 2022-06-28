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

"""
Calculate interactions using ODDT (https://oddt.readthedocs.io/).

Example usage:
    $ python -m pipelines.xchem.calc_interactions -p ../../data/mpro/Mpro-x0387_0.pdb -i ../../data/mpro/hits-17.sdf.gz -o output

The interactions that are searched for are:
* Hydrogen bond
* Halogen bond
* Hydrophobic
* Salt bridge
* Pi stacking
* Pi cation
See the ODDT docs for more details.

We introduce on additional concept, the 'canonical' interaction at the protein.
This allows the same interaction to be identified when comparing structures of the same protein.
The canonical representation comprise the protein residue sometimes followed in the case of H-bond interactions by a
definition of which part of that residue is involved (BN: backbone N, BO: backbone O, SC, sidechain). Other types of
interaction do not need this distinction so just use the residue type and number.

In the case of hydrophobic interactions there can be multiple interactions between a ligand and a residue. The information
recorded is only for the shortest of these interactions.

The different interaction types are written as data fields with a syntax like:

GLN189SC [10.922, 4.777, 26.595] [12.074, 5.388, 23.622] 3.246 3
HIS41SC [12.815, -4.64, 20.946] [13.769, -1.898, 21.404] 2.939 15

The tokens are:
1. The canonical site on the protein
2. The 3D location of the protein atom or pseudo atom as x, y, z coordinates inside square brackets
3. The 3D location of the ligand atom or pseudo atom as x, y, z coordinates inside square brackets
4. The distance between the atoms or pseudo atoms
5. Atom number (zero based) of the ligand atom. This is missing where the ligand part is a pseudo atom.

Normally the interaction is between 2 atoms, but in some cases like Pi interactions it is to an aromatic ring that is
represented as a pseudo atom at the centre of the ring.

You can also specify 'key interactions' that are counted. e.g in the case of the above H-bond interactions you can use
the `--key-hbond GLN189SC HIS41SC' option to count those specific interactions. One or more values are allowed. Similar
things for other types of interaction using the `--key-*` options.

Typically you specify a single protein and a set of ligands (e.g. as a SDF) to work on. Each ligand is calculated against
that protein. However, you can also specify multiple proteins as the --protein option allows multiple values.
In this case the proteins are used in turn (first ligand uses the first protein, second ligand uses the second protein ...)
until the proteins run out after which the last one continues to be used. In the usual case just specify a single protein
and it will be used for all ligands.

The following data fields will be generated. Not all will be present depending of the options used and whether particular
interactions are found:
* HydrophobicInteractions - details of the H-bond interactions
* HalogenBondInteractions - details of the halogen bond interactions
* HydrophobicInteractions - details of the hydrophobic interactions
* SaltBridgeInteractions - details of the salt bridge interactions
* PiStackingInteractions - details of the Pi stacking interactions
* PiCationInteractions - details of the Pi cation interactions
* NumTotalInteractions - total number of interactions found
* NumKeyInteractions - number of key interactions found
* KeyInteractions - the key interactions that were found
* vina_* - a whole series of AutoDock VINA scores that are generated when running NNScore
"""

from __future__ import print_function
import argparse, traceback
import json, sys

import utils, interact
from dm_job_utilities.dm_log import DmLog

import oddt
from oddt import toolkit, spatial, interactions

from rdkit import Chem


# start function definitions #########################################

def generate_js(type, values):
    content = []
    for v in values:
        js = "  ['%s', '%s', [%s, %s, %s], [%s, %s, %s]]," % (
            type, v[0], v[1][0], v[1][1], v[1][2], v[2][0], v[2][1], v[2][2])
        content.append(js)
    return "\n".join(content)


def read_next_protein(proteins, format, previous, index, keep_hs=False):

    if previous and index >= len(proteins):
        return previous
    protein = next(toolkit.readfile(format, proteins[index]))
    if not protein:
        raise ValueError('Unable to read protein')
    else:
        utils.log('Read protein', index + 1)
        if not keep_hs:
            protein.removeh()
        protein.protein = True
        return protein


def process(protein_files, ligands, output_file, key_inters, filter_strict=False,
            exact_protein=False, exact_ligand=False, keep_hs_protein=False, keep_hs_ligand=False, report_file=None,
            compare_file=None, interval=1000):

    if len(protein_files) > 1:
        utils.log(len(protein_files), 'proteins specified')

    utils.expand_path(output_file)

    if report_file:
        report_data = []
        utils.expand_path(report_file)
    else:
        report_data = None

    if compare_file:
        with open(compare_file, "r") as f:
            txt = f.read()
            compare_data = interact.from_json(txt)
    else:
        compare_data = None

    total = 0
    count = 0
    errors = 0
    protein = None

    with Chem.SDWriter(output_file) as writer:

        for ligand in oddt.toolkit.readfile('sdf', ligands):
            # utils.log('Processing ligand', total + 1)
            try:
                if not keep_hs_ligand:
                    ligand.removeh()
                protein = read_next_protein(protein_files, 'pdb', protein, total, keep_hs=keep_hs_protein)
                inter_data = interact.process_mol(protein, ligand, key_inters, total, filter_strict=filter_strict,
                                         exact_protein=exact_protein, exact_ligand=exact_ligand, compare_data=compare_data)
                if report_data is not None:
                    report_data.append(inter_data)

                # write the RDKit mol
                writer.write(ligand.Mol)
                count += 1
            except KeyboardInterrupt:
                utils.log('Interrupted')
                sys.exit(0)
            except:
                errors += 1
                traceback.print_exc()
                sys.exit(0)
            finally:
                total += 1

            if interval and total % interval == 0:
                DmLog.emit_event("Processed {} molecules".format(total))

        # print(json.dumps(report_data, cls=interact.InteractionEncoder))

        if report_data:
            with open(report_file, 'w') as report:
                json.dump(report_data, report, cls=interact.InteractionEncoder)

    return count, errors


### start main execution #########################################

def main():
    # Example usage
    # python oddt_interactions.py -p data/dhfr-receptor.pdb -i data/dhfr-ligand.mol -o output.sdf

    parser = argparse.ArgumentParser(description='Calculate interactions')
    parser.add_argument('-i', '--input', required=True, help="SD file with poses to analyse")
    parser.add_argument('-p', '--protein', nargs='+', required=True, help="File with protein (PDB conda activate format")
    parser.add_argument('-o', '--output', required=True, help="SD file for output")
    parser.add_argument('--strict', action='store_true', help='Strict filtering')
    parser.add_argument('--exact-protein', action='store_true', help='Exact matching of hydrogens and charges for protein')
    parser.add_argument('--exact-ligand', action='store_true', help='Exact matching of hydrogens and charges for ligand')
    parser.add_argument('--keep-hs-protein', action='store_true', help='Keep hydrogens on the protein')
    parser.add_argument('--keep-hs-ligand', action='store_true', help='Keep hydrogens on the ligand')
    parser.add_argument('--key-hbond', nargs='*', help='List of canonical H-bond interactions to count')
    parser.add_argument('--key-hydrophobic', nargs='*', help='List of canonical hydrophobic interactions to count')
    parser.add_argument('--key-salt-bridge', nargs='*', help='List of canonical salt bridge interactions to count')
    parser.add_argument('--key-pi-stacking', nargs='*', help='List of canonical pi stacking interactions to count')
    parser.add_argument('--key-pi-cation', nargs='*', help='List of canonical pi cation interactions to count')
    parser.add_argument('--key-halogen', nargs='*', help='List of canonical halogen bond interactions to count')
    parser.add_argument('-r', '--report-file', help="File for the report")
    parser.add_argument('-c', '--compare', help="Compare interactions with this report")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log("Calculate interactions Args: ", args)

    key_inters = {}
    if args.key_hbond:
        key_inters[interact.I_TYPE_HBOND] = args.key_hbond
    if args.key_hydrophobic:
        key_inters[interact.I_TYPE_HYDROPHOBIC] = args.key_hydrophobic
    if args.key_salt_bridge:
        key_inters[interact.I_TYPE_SALT_BRIDGE] = args.key_salt_bridge
    if args.key_pi_stacking:
        key_inters[interact.I_TYPE_PI_STACKING] = args.key_pi_stacking
    if args.key_pi_cation:
        key_inters[interact.I_TYPE_PI_CATION] = args.key_pi_cation
    if args.key_halogen:
        key_inters[interact.I_TYPE_HALOGEN] = args.key_halogen


    # this does the processing
    count, errors = process(args.protein, args.input, args.output, key_inters, filter_strict=args.strict,
                            exact_protein=args.exact_protein, exact_ligand=args.exact_ligand,
                            keep_hs_protein=args.keep_hs_protein, keep_hs_ligand=args.keep_hs_ligand,
                            report_file=args.report_file, compare_file=args.compare, interval=args.interval)
    DmLog.emit_event('Processing complete.', count, 'records processed.', errors, 'errors')
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
