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
Reads a file of SMILES and generates low energy conformers for each molecule using the "sharded file system".

See also le_conformers_for_mol.py that provides a simpler way to generate conformers for a single molecule using the
same methodology that this module uses.
"""
import os, argparse, traceback, time, gzip
import utils

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors


def determine_num_confs(mol):
    """
    Determine the number of conformers to generate.

    These rules are defined in Ejeber et al.J. Chem. Inf. Model. 2012, 52, 1146-1158
    https://pubs.acs.org/doi/abs/10.1021/ci2004658
    Note that we use a different RMS threshold as the default, 1.0 as opposed to 0.35

    :param mol: The molecule
    :return: The number of conformers to generate
    """
    rotb = Descriptors.NumRotatableBonds(mol)
    if rotb <= 7:
        return 50
    elif rotb <= 12:
        return 200
    else:
        return 300


def gen_conformers(mol, rms_threshold, minimize_cycles, remove_hydrogens):

    num_confs = determine_num_confs(mol)
    molh = Chem.AddHs(mol)

    cids = AllChem.EmbedMultipleConfs(molh, numConfs=num_confs)

    res = AllChem.MMFFOptimizeMoleculeConfs(molh, maxIters=minimize_cycles)
    # energies_with_index will contain tuples of (conf_index, not_converged, energy)
    energies_with_index = []
    for idx, val in enumerate(res):
        energies_with_index.append((idx, val[0], val[1]))
    # sort by energy
    energies_with_index.sort(key=lambda t: t[2])

    # set properties on the conformers
    base_energy = energies_with_index[0][2]
    for val in energies_with_index:
        energy = val[2]
        energy_delta = energy - base_energy
        conf = molh.GetConformer(val[0])
        conf.SetDoubleProp('Energy', energy)
        conf.SetDoubleProp('Energy_Delta', energy_delta)

    to_keep = [energies_with_index[0]]
    utils.log('  Added', to_keep[0])
    for i in range(1, len(energies_with_index)):
        to_test = energies_with_index[i]
        #utils.log('Checking', i, to_test)
        lowest_rms = 999999999
        for val in to_keep:
            rms = AllChem.GetConformerRMS(molh, to_test[0], val[0])
            if rms < lowest_rms:
                lowest_rms = rms
        if lowest_rms > rms_threshold:
            to_keep.append(to_test)
        #     utils.log('  Added', to_test, lowest_rms)
        # else:
        #     utils.log('  Skipping', to_test, lowest_rms)

    if remove_hydrogens:
        molh = Chem.RemoveHs(molh)

    # create new mol for the conformers we want to keep
    final_mol = Chem.RWMol(molh)
    final_mol.RemoveAllConformers()

    for item in to_keep:
        idx = item[0]
        final_mol.AddConformer(molh.GetConformer(idx), assignId=True)

    utils.log("Num conformers:", molh.GetNumConformers(), '->', final_mol.GetNumConformers())
    return final_mol


def execute(input, data_dir, minimize_cycles=500, remove_hydrogens=False, rms_threshold=1.0, interval=None):

    utils.log_dm_event('Executing ...')

    input_count = 0
    enumerated_count = 0
    conformer_count = 0
    error_count = 0

    with open(input) as infile:
        for line in infile:
            input_count += 1
            conf_count_for_mol = 0
            try:
                tokens = line.strip().split('\t')
                std_smi = tokens[0]
                uid = tokens[1]
                digest = tokens[2]

                if interval and input_count % interval == 0:
                    utils.log_dm_event("Processed {} records".format(input_count))

                parts = [data_dir]
                parts.extend(utils.get_path_from_digest(digest))
                path = os.path.join(*parts)
                if not os.path.isdir(path):
                    utils.log_dm_event('WARNING, path', path, 'not found')
                    error_count += 1
                    continue

                smi_in = os.path.join(path, digest + '.smi')
                sdf_out = os.path.join(path, digest + '_le_confs.sdf.gz')
                if not os.path.exists(smi_in):
                    utils.log_dm_event('WARNING, smiles file', smi_in, 'not found')
                    error_count += 1
                    continue

                utils.log('INFO, Processing file', smi_in)
                with gzip.open(sdf_out, 'wt') as gz:
                    with Chem.SDWriter(gz) as writer:
                        with open(smi_in) as enums:
                            for line in enums:
                                enumerated_count += 1
                                tokens2 = line.strip().split('\t')
                                enum_smi = tokens2[0]
                                uid2 = tokens2[1]
                                code = tokens2[2]

                                mol = Chem.MolFromSmiles(enum_smi)
                                if not mol:
                                    error_count += 1
                                    utils.log_dm_event("ERROR, Failed to create molecule", input_count, enumerated_count)
                                    continue

                                molh = Chem.AddHs(mol)
                                molh.SetProp('_Name', uid2)
                                molh.SetProp('std_smi', std_smi)
                                molh.SetProp('enum_smi', enum_smi)
                                molh.SetProp('enum_code', code)
                                if code != 'B':
                                    molh.SetProp('parent_uuid', uid)

                                mol_with_confs = gen_conformers(molh, rms_threshold, minimize_cycles, remove_hydrogens)
                                for idx in range(mol_with_confs.GetNumConformers()):
                                    mol_with_confs.SetDoubleProp('Energy', mol_with_confs.GetConformer(idx).GetDoubleProp('Energy'))
                                    mol_with_confs.SetDoubleProp('Energy_Delta', mol_with_confs.GetConformer(idx).GetDoubleProp('Energy_Delta'))
                                    writer.write(mol_with_confs, confId=idx)
                                    conformer_count += 1
                                    conf_count_for_mol += 1
                                mol_with_confs.ClearProp('Energy')
                                mol_with_confs.ClearProp('Energy_Delta')

                utils.log('INFO, Generated', conf_count_for_mol, 'conformers for', std_smi)

            except:
                error_count += 1
                traceback.print_exc()

    return input_count, enumerated_count, conformer_count, error_count


### start main execution #########################################

def main():

    # Example:
    #   python3 le_conformers.py -i bar.smi

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate conformers')
    parser.add_argument('-i', '--input', required=True, help="Input file as SMILES")
    parser.add_argument('--data-dir', default='molecules/sha256', help="Data directory")
    parser.add_argument('-m', '--minimize-cycles', type=int, default=500, help="Number of MMFF minimisation cycles")
    parser.add_argument('-t', '--rms-threshold', type=float, default=1.0, help="RMS threshold for excluding conformers")
    parser.add_argument('--remove-hydrogens', action='store_true', help='Remove hydrogens from the output')
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log_dm_event("le_conformers: ", args)

    start = time.time()
    input_count, enumerated_count, conformer_count, error_count = \
        execute(args.input, args.data_dir,
                minimize_cycles=args.minimize_cycles,
                remove_hydrogens=args.remove_hydrogens,
                rms_threshold=args.rms_threshold,
                interval=args.interval
                )
    end = time.time()

    utils.log_dm_event('Inputs:', input_count, 'Enumerated:', enumerated_count,
                       'Conformers:', conformer_count, 'Errors:', error_count, 'Time (s):', end - start)


if __name__ == "__main__":
    main()
