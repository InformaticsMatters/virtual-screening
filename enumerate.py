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
Reads a file of SMILES and enumerates undefined chiral centres, tautomers and charge states.
Writes the enumerated molecules as SMILES and a single 3D conformer as SDF.

Note that RDKit does not write the charge information to the atom block (this is deprecated in favour of using
'M   CHG' records), but some old software like rDock requires the charges to be present in the atom block. This
script patches the atom blocks to add the charge information. See rdkit_utils.py for details.


"""
import os, sys, argparse, traceback, uuid, gzip
import utils, rdkit_utils
from dm_job_utilities.dm_log import DmLog

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

from dimorphite_dl import run_with_mol_list


def enumerate_undefined_chirals(mol, tryEmbedding=False):
    """
    Enumerate undefined chiral centres that are points of substitution.
    Chiral centres that are already specifically defined are left untouched.
    :param mol: The mol to enumerate
    :return: A list of enumerated molecules
    """
    opts = StereoEnumerationOptions(tryEmbedding=tryEmbedding)
    isomers = tuple(EnumerateStereoisomers(mol, options=opts))
    return isomers


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
                if atom.GetFormalCharge() != 0:
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
        mol.SetProp('enum_code', code)
        mydict[smiles] = mol


def execute(input, data_dir, delimiter='\t',
            enumerate_chirals=False,
            enumerate_charges=False, enumerate_tautomers=False,
            combinatorial=False,
            min_charge=None, max_charge=None, num_charges=None,
            min_hac=None, max_hac=None,
            min_ph=5.0, max_ph=9.0,
            add_hydrogens=True,
            try_embedding=False,
            max_tautomers=25,
            interval=0):

    DmLog.emit_event('Executing ...')

    if enumerate_tautomers or combinatorial:
        enumerator = TautomerEnumerator()
        enumerator.SetMaxTautomers(max_tautomers)

    count = 0
    total = 0
    errors = 0
    excluded = 0
    
    with open(input) as infile:
        for line in infile:
            count += 1
            
            tokens = line.strip().split(delimiter)
            smi = tokens[0]
            uid = tokens[1]
            digest = tokens[2]
            
            if interval and count % interval == 0:
                    DmLog.emit_event("Processed {} records".format(count))
                    
            mol = Chem.MolFromSmiles(smi)
            if not mol:
                errors += 1
                DmLog.emit_event("Failed to create molecule", count)
                continue

            # if we have HAC filters apply them
            hac = mol.GetNumHeavyAtoms()
            if min_hac is not None and min_hac > hac:
                excluded += 1
                continue
            if max_hac is not None and max_hac < hac:
                excluded += 1
                continue
            
            parts = [data_dir]
            parts.extend(utils.get_path_from_digest(digest))
            path = os.path.join(*parts, digest)
            #utils.log('Handling', path)

            with open(path + '.smi', 'w') as smi_writer:
                with gzip.open(path + '.sdf.gz', 'wt') as gz:

                    try:
                        # use a dict as we want to keep the order
                        enumerated_mols = {}
                        add_molecule(enumerated_mols, mol, 'B')

                        if enumerate_chirals:
                            mols = enumerate_undefined_chirals(mol, tryEmbedding=try_embedding)
                            add_molecules(enumerated_mols, mols, 'C')

                        if combinatorial:
                            tautomers = gen_tautomers(enumerated_mols.values(), enumerator)
                            protonated_mols = gen_charges(tautomers, min_ph, max_ph, min_charge, max_charge,
                                                          num_charges)
                            add_molecules(enumerated_mols, protonated_mols, 'X')
                        else:
                            if enumerate_tautomers:
                                tautomers = gen_tautomers(enumerated_mols.values(), enumerator)
                            if enumerate_charges:
                                protonated_mols = gen_charges(enumerated_mols.values(), min_ph, max_ph, min_charge,
                                                              max_charge, num_charges)
                            if enumerate_tautomers:
                                add_molecules(enumerated_mols, tautomers, 'T')
                            if enumerate_charges:
                                add_molecules(enumerated_mols, protonated_mols, 'M')

                        total += len(enumerated_mols)

                        for m in enumerated_mols.values():

                            code = m.GetProp('enum_code')
                            if code == 'B':
                                u = uid
                            else:
                                u = str(uuid.uuid4())

                            enum_smi = Chem.MolToSmiles(m)

                            # generate 3D conformer
                            m2 = Chem.AddHs(m)
                            AllChem.EmbedMolecule(m2)
                            AllChem.MMFFOptimizeMolecule(m2)
                            if not add_hydrogens:
                                m2 = Chem.RemoveHs(m2)

                            m2.SetProp('uuid', u)
                            m2.SetProp('_Name', u)
                            m2.SetProp('sha256', digest)
                            m2.SetProp('std_smi', smi)
                            m2.SetProp('enum_smi', enum_smi)

                            # write to SDF
                            sdf_block = Chem.SDWriter.GetText(m2)
                            chg_block = rdkit_utils.updateChargeFlagInAtomBlock(sdf_block)
                            gz.write(chg_block)

                            # write the SMILES last so that it's only written if the 3D generation is successful
                            smi_writer.write(enum_smi + '\t' + u + '\t' + code + '\n')

                    except:
                        errors += 1
                        traceback.print_exc()


    return count, total, excluded, errors


### start main execution #########################################

def main():

    # Example:
    #   ./enumerate.py -i bar.smi --enumerate-tautomers --enumerate-chirals --enumerate-charges

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate candidates')
    parser.add_argument('-i', '--input', required=True, help="Input file as SMILES")
    parser.add_argument('--data-dir', default='molecules/sha256', help="Data directory")
    parser.add_argument('-d', '--delimiter', default='\t', help="Delimiter")
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
    parser.add_argument('--try-embedding', action='store_true', help='Try to embed a 3D molecule to verify the stereochemisty is sane')
    parser.add_argument("--max-tautomers", type=int, default=25, help="Maximum number of tautomers to generate")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("enumerate_candidates: ", args)

    # save the arguments
    input = args.input
    data_dir = args.data_dir
    delimiter = args.delimiter
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
    try_embedding = args.try_embedding
    max_tautomers = args.max_tautomers
    interval = args.interval 


    # Dimporphite needs to use argparse with its own arguments, not messed up with our arguments
    # so we store the original args
    orig_sys_argv = sys.argv[:]

    # Remove all the parameters, keeping only the filename (first one) so that
    # dimorphite is unaware of any previous commandline parameters.
    sys.argv = sys.argv[:1]

    count, total, excluded, errors = execute(input, data_dir, delimiter=delimiter,
                                             enumerate_chirals=enumerate_chirals,
                                             enumerate_charges=enumerate_charges,
                                             enumerate_tautomers=enumerate_tautomers,
                                             combinatorial=combinatorial,
                                             min_charge=min_charge, max_charge=max_charge, num_charges=num_charges,
                                             min_hac=min_hac, max_hac=max_hac,
                                             min_ph=min_ph, max_ph=max_ph,
                                             add_hydrogens=add_hydrogens,
                                             try_embedding=try_embedding,
                                             max_tautomers=max_tautomers,
                                             interval=interval
                                             )

    DmLog.emit_event('Count:', count, 'Total', total, 'Excluded:', excluded, 'Errors:', errors)



if __name__ == "__main__":
    main()
