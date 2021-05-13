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


"""
Reads a file of SMILES and enumerates undefined chiral centres, tautomers and charge states.
Writes results as SMILES.
"""
import os, sys, argparse, traceback
import utils

from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator

from dimorphite_dl import run_with_mol_list


def enumerate_undefined_chirals(mol):
    """
    Enumerate undefined chiral centres that are points of substitution.
    Chiral centres that are already specifically defined are left untouched.
    :param mol: The mol to enumerate
    :return:
    """

    global num_swap_success

    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    #utils.log('  Chiral centers:', chiral_centers, 'Free atom counts:', free_counts)

    chiral_targets = []
    for i, cc in enumerate(chiral_centers):
        if cc[1] == '?':  # only enumerate if the stereo is undefined
            chiral_targets.append(cc[0])

    if chiral_targets:
        chiral_swaps = []
        mol_cw = Chem.RWMol(mol)
        mol_cw.GetAtomWithIdx(chiral_targets[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        chiral_swaps.append(mol_cw)
        mol_ccw = Chem.RWMol(mol)
        mol_ccw.GetAtomWithIdx(chiral_targets[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        chiral_swaps.append(mol_ccw)
        for j in range(1, len(chiral_targets)):
            new_enum_mols = []
            for m in chiral_swaps:
                nm1 = Chem.RWMol(m)
                nm1.GetAtomWithIdx(chiral_targets[j]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                new_enum_mols.append(nm1)
                nm2 = Chem.RWMol(m)
                nm2.GetAtomWithIdx(chiral_targets[j]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                new_enum_mols.append(nm2)
                chiral_swaps = new_enum_mols
                num_swap_success += len(chiral_swaps)
    else:
        chiral_swaps = [Chem.RWMol(mol)]
    #utils.log('  Mols to embed:', len(chiral_swaps))
    return chiral_swaps

    
def gen_tautomers(mols, enumerator):
    tautomers = []
    for m in mols:
        ts = enumerator.Enumerate(m)
        tautomers.extend(ts)
    return tautomers

def gen_charges(mols, min_ph, max_ph, min_charge, max_charge, num_charges):
    protonated_mols = run_with_mol_list(mols, min_ph=min_ph, max_ph=max_ph)
    wanted_mols = []
    for mol in protonated_mols:
        if num_charges is not None:
            nc = 0
            for atom in mol.GetAtoms():
                if atom.GetFormalCharge() > 0:
                    nc += 1
            if nc > num_charges:
                continue

        if min_charge is not None or max_charge is not None:
            charge = Chem.GetFormalCharge(mol)
            if min_charge is not None and charge < min_charge:
                continue
            if max_charge is not None and charge > max_charge:
                continue
        wanted_mols.append(mol)
    return wanted_mols

def add_molecules(mydict, mols, code):
    for mol in mols:
        add_molecule(mydict, mol, code)


def add_molecule(mydict, mol, code):
    smiles = Chem.MolToSmiles(mol)
    if smiles not in mydict:
        mol.SetProp('CODE', code)
        mydict[smiles] = mol


def execute(suppl, data_dir, enumerate_chirals=False,
            enumerate_charges=False, enumerate_tautomers=False,
            combinatorial=False,
            min_charge=None, max_charge=None, num_charges=None,
            min_hac=None, max_hac=None,
            min_ph=5.0, max_ph=9.0,
            add_hydrogens=False,
            interval=0):

    utils.log('Executing ...')

    if enumerate_tautomers or combinatorial:
        enumerator = TautomerEnumerator()

    count = 0
    total = 0
    errors = 0
    excluded = 0
    for mol in suppl:
        count += 1
        
        if interval and count % interval == 0:
                utils.log("Processed {} records".format(count))

        if not mol:
            errors += 1
            utils.log("Failed to create molecule", count)
            continue

        # if we have HAC filters apply them
        hac = mol.GetNumHeavyAtoms()
        if min_hac is not None and min_hac > hac:
            excluded += 1
            continue
        if max_hac is not None and max_hac < hac:
            excluded += 1
            continue

        digest = mol.GetProp('_Name')
        #utils.log('Processing mol', count)
        
        parts = [data_dir]
        parts.extend(utils.get_path_from_digest(digest))
        path = os.path.join(*parts, digest + '.smi')
        #utils.log('Handling', path)
        
        codes = ''
        with open(path, 'w') as writer:

            try:

                # use a dict as we want to keep the order
                enumerated_mols = {}
                add_molecule(enumerated_mols, mol, 'B')

                if enumerate_chirals:
                    mols = enumerate_undefined_chirals(mol)
                    add_molecules(enumerated_mols, mols, 'C')

                if combinatorial:
                    tautomers = gen_tautomers(enumerated_mols.values(), enumerator)
                    protonated_mols = gen_charges(tautomers, min_ph, max_ph, min_charge, max_charge, num_charges)
                    add_molecules(enumerated_mols, protonated_mols, 'X')
                else:
                    if enumerate_tautomers:
                        tautomers = gen_tautomers(enumerated_mols.values(), enumerator)
                    if enumerate_charges:
                        protonated_mols = gen_charges(enumerated_mols.values(), min_ph, max_ph, min_charge, max_charge, num_charges)
                    if enumerate_tautomers:
                        add_molecules(enumerated_mols, tautomers, 'T')
                    if enumerate_charges:
                        add_molecules(enumerated_mols, protonated_mols, 'M')

                total += len(enumerated_mols)

                for m in enumerated_mols.values():
                    writer.write(Chem.MolToSmiles(m) + '\t' + m.GetProp('CODE') + '\n')

            except:
                errors += 1
                traceback.print_exc()


    return count, total, excluded, errors


### start main execution #########################################

def main():

    # Example:
    #   python3 enumerate.py -i bar.smi --data-dir combined --enumerate-tautomers --enumerate-chirals --enumerate-charges

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate candidates')
    parser.add_argument('-i', '--input', required=True, help="Input file as SMILES")
    parser.add_argument('--data-dir', required=True, help="Data directory")
    parser.add_argument('-d', '--delimiter', default='\t', help="Delimiter")
    parser.add_argument("--name-column", type=int, default=1, help="Column for name field (zero based)")
    parser.add_argument('-t', '--title-line', default=False, help="Do the files have title lines")
    parser.add_argument('--enumerate-charges', help='Enumerate charge forms', action='store_true')
    parser.add_argument('--enumerate-chirals', help='Enumerate undefined chiral centers', action='store_true')
    parser.add_argument('--enumerate-tautomers', help='Enumerate undefined chiral centers', action='store_true')
    parser.add_argument('--combinatorial', help='Combinatorial enumeration of charge and tautomers', action='store_true')
    parser.add_argument("--min-hac", type=int, help="Minimum heavy atom count to consider")
    parser.add_argument("--max-hac", type=int, help="Maximum heavy atom count to consider")
    parser.add_argument("--min-ph", type=float, default=5, help="Minimum pH for charge enumeration")
    parser.add_argument("--max-ph", type=float, default=9, help="Maximum pH for charge enumeration")
    parser.add_argument("--min-charge", type=int, help="Minimum charge of molecule to process")
    parser.add_argument("--max-charge", type=int, help="Maximum charge of molecule to process")
    parser.add_argument("--num-charges", type=int, help="Maximum number of atoms with a charge")
    parser.add_argument('--add-hydrogens', action='store_true', help='Include hydrogens in the output')
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log("enumerate_candidates: ", args)

    # save the arguments
    input = args.input
    data_dir = args.data_dir
    delimiter = args.delimiter
    title_line = args.title_line
    name_column = args.name_column
    enumerate_charges = args.enumerate_charges
    enumerate_chirals = args.enumerate_chirals
    enumerate_tautomers = args.enumerate_tautomers
    combinatorial = args.combinatorial
    min_hac = args.min_hac
    max_hac = args.max_hac
    min_ph = args.min_ph
    max_ph = args.max_ph
    min_charge = args.min_charge
    max_charge = args.max_charge
    num_charges = args.num_charges
    add_hydrogens = args.add_hydrogens
    interval = args.interval 


    # Dimporphite needs to use argparse with its own arguments, not messed up with our arguments
    # so we store the original args
    orig_sys_argv = sys.argv[:]

    # Remove all the parameters, keeping only the filename (first one) so that
    # dimorphite is unaware of any previous commandline parameters.
    sys.argv = sys.argv[:1]

    suppl = Chem.SmilesMolSupplier(input, delimiter=delimiter, titleLine=title_line, nameColumn=name_column)

    count, total, excluded, errors = execute(suppl, data_dir, enumerate_chirals=enumerate_chirals,
                                   enumerate_charges=enumerate_charges, enumerate_tautomers=enumerate_tautomers,
                                   combinatorial=combinatorial,
                                   min_charge=min_charge, max_charge=max_charge, num_charges=num_charges,
                                   min_hac=min_hac, max_hac=max_hac,
                                   min_ph=min_ph, max_ph=max_ph, add_hydrogens=add_hydrogens,
                                   interval=interval
                                   )

    utils.log(count, total, excluded, errors)



if __name__ == "__main__":
    main()
