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

import argparse, time, traceback
from dm_job_utilities.dm_log import DmLog

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdFMCS

import utils


def align(query_mol, inputs_sdf, outputs_sdf, query_atoms=None, keep_all=False, mcs_params={}, interval=None):
    count = 0
    errors = 0
    if query_mol.endswith('.mol'):
        query = Chem.MolFromMolFile(query_mol)
    elif query_mol.endswith('.sdf'):
        # read the first mol
        with Chem.SDMolSupplier(query_mol) as supplr:
            query = next(supplr)
    else:
        DmLog.emit_event("Unsupported query format. Must be .mol or .sdf. Found:", query_mol)
        exit(1)

    if not query:
        DmLog.emit_event("Failed to read query molecule. Can't continue")
        exit(1)

    utils.log('Using MCS params:', mcs_params)

    with Chem.SDMolSupplier(inputs_sdf) as supplr:
        with Chem.SDWriter(outputs_sdf) as writer:
            for mol in supplr:
                count += 1
                if not mol:
                    errors += 1
                    utils.log("Failed to read molecule", count)
                    continue
                try:
                    smartsmol = find_mcs([mol, query], **mcs_params)
                    mappings = find_mappings(mol, query, smartsmol, query_atoms)
                    best_rmsd = 999999999999
                    best_mol = None
                    for mapping in mappings:
                        rmsd = align_mol(mol, query, atom_mapping=mapping)
                        if rmsd < best_rmsd:
                            best_rmsd = rmsd
                            best_mol = mol
                        if keep_all:
                            writer.write(mol)
                    if best_mol and not keep_all:
                        writer.write(mol)
                except:
                    errors += 1
                    DmLog.emit_event("Error processing molecule", count)
                    traceback.print_exc()

                if interval and count % interval == 0:
                    DmLog.emit_event("Processed {} records, {} errors".format(count, errors))

    return count, errors


def align_mol(probe, reference, atom_mapping=[]):
    rmsd = rdMolAlign.AlignMol(probe, reference, atomMap=atom_mapping)
    probe.SetDoubleProp('RMSD', rmsd)
    return rmsd


def find_mappings(probe, reference, smartsmol, query_atoms):
    matches1 = reference.GetSubstructMatches(smartsmol, uniquify=True)
    matches2 = probe.GetSubstructMatches(smartsmol, uniquify=True)
    utils.log("Num matches:", len(matches1), len(matches2))
    utils.log(matches1)
    utils.log(matches2)
    mappings = []
    for match1 in matches1:
        for match2 in matches2:
            mapping = []
            mappings.append(mapping)
            for a, b in zip(match1, match2):
                if not query_atoms or a in query_atoms:
                    mapping.append((b, a))
    utils.log(mappings)
    return mappings


def find_mcs(mols, maximizeBonds=True, matchValences=False, ringMatchesRingOnly=False,
             completeRingsOnly=False, matchChiralTag=False, atomCompare='CompareElements',
             bondCompare='CompareOrder', ringCompare='IgnoreRingFusion'):

    atomCompareVal = rdFMCS.AtomCompare.names[atomCompare]
    bondCompareVal = rdFMCS.BondCompare.names[bondCompare]
    ringCompareVal = rdFMCS.RingCompare.names[ringCompare]
    mcs = rdFMCS.FindMCS(mols, maximizeBonds=maximizeBonds, matchValences=matchValences,
                         ringMatchesRingOnly=ringMatchesRingOnly, completeRingsOnly=completeRingsOnly,
                         matchChiralTag=matchChiralTag, atomCompare=atomCompareVal, bondCompare=bondCompareVal,
                         ringCompare=ringCompareVal)
    smartsmol = Chem.MolFromSmarts(mcs.smartsString)
    return smartsmol


### start main execution #########################################

def main():

    # Example:
    #   ./align_mol.py -i conformers.sdf -q query.mol -o aligned.sdf

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate conformers for molecule')
    parser.add_argument('-i', '--input', required=True, help="Input molecules as SDF")
    parser.add_argument('-q', '--query', required=True, help="Query molecule as molfile")
    parser.add_argument('-o', '--output', required=True, help="Output as SD file")
    parser.add_argument('-a', '--align-atoms',
                        help="Align using only these atoms from query (comma separated list, zero based)")
    parser.add_argument('-k', '--keep-all', action='store_true', help='Keep all the matches rather than just the best')
    # For details of the MCS parameters look at:
    # http://rdkit.org/docs/source/rdkit.Chem.rdFMCS.html?highlight=findmcs#rdkit.Chem.rdFMCS.FindMCS
    parser.add_argument('--mcs-maximize-bonds', type=bool, help="maximizeBonds param for MCS")
    parser.add_argument('--mcs-match-valences', type=bool, help="matchValences param for MCS")
    parser.add_argument('--mcs-ring-matches-ring-only', type=bool, help="ringMatchesRingOnly param for MCS")
    parser.add_argument('--mcs-complete-rings-only', type=bool, help="completeRingsOnly param for MCS")
    parser.add_argument('--mcs-match-chiral-tag', type=bool, help="matchChiralTag param for MCS")
    parser.add_argument('--mcs-atom-compare', help="bondCompare param for MCS")
    parser.add_argument('--mcs-bond-compare', help="atomCompare param for MCS")
    parser.add_argument('--mcs-ring-compare', help="ringCompare param for MCS")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("align_mol: ", args)

    if args.align_atoms:
        query_atoms = list(map(lambda s: int(s), args.align_atoms.split(',')))
        utils.log('Query atoms:', query_atoms)
    else:
        query_atoms = None

    mcs_params = {}
    if args.mcs_maximize_bonds:
        mcs_params['maximizeBonds'] = args.mcs_maximize_bonds
    if args.mcs_match_valences:
        mcs_params['matchValences'] = args.mcs_match_valences
    if args.mcs_ring_matches_ring_only:
        mcs_params['ringMatchesRingOnly'] = args.mcs_ring_matches_ring_only
    if args.mcs_complete_rings_only:
        mcs_params['completeRingsOnly'] = args.mcs_complete_rings_only
    if args.mcs_match_chiral_tag:
        mcs_params['matchChiralTag'] = args.mcs_match_chiral_tag
    if args.mcs_atom_compare:
        mcs_params['atomCompare'] = args.mcs_atom_compare
    if args.mcs_bond_compare:
        mcs_params['bondCompare'] = args.mcs_bond_compare
    if args.mcs_ring_compare:
        mcs_params['ringCompare'] = args.mcs_ring_compare

    start = time.time()
    count, errors = align(args.query, args.input, args.output,
                          query_atoms=query_atoms, keep_all=args.keep_all, mcs_params=mcs_params, interval=args.interval)
    end = time.time()

    DmLog.emit_event('Processed', count, 'molecules.', errors, 'errors. Time (s):', end - start)
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
