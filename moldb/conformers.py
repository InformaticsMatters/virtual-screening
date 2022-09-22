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
Reads a file of SMILES and generates low energy conformers for each molecule using the "sharded file system".

See also le_conformers_for_mol.py that provides a simpler way to generate conformers for a single molecule using the
same methodology that this module uses.
"""
import os, argparse, traceback, time, gzip
import utils
from dm_job_utilities.dm_log import DmLog

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


def gen_conformers(mol, rms_threshold, minimize_cycles, remove_hydrogens, num_conformers=None):

    if num_conformers:
        num_confs = num_conformers
    else:
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
    #utils.log('  Added', to_keep[0])
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

    #utils.log("Num conformers:", molh.GetNumConformers(), '->', final_mol.GetNumConformers())
    return final_mol


def execute(input, output, minimize_cycles=500, rms_threshold=1.0, num_conformers=None, interval=None):

    DmLog.emit_event('Executing ...')

    input_count = 0
    enumerated_count = 0
    conformer_count = 0
    error_count = 0

    utils.expand_path(output)

    with open(input) as infile:
        with open(output, 'wt') as writer:
            for line in infile:
                input_count += 1
                conf_count_for_mol = 0
                try:
                    tokens = line.strip().split('\t')
                    smi = tokens[0]
                    id = tokens[1]

                    if interval and input_count % interval == 0:
                        DmLog.emit_event("Processed {} records".format(input_count))

                    mol = Chem.MolFromSmiles(smi)
                    if not mol:
                        error_count += 1
                        DmLog.emit_event("ERROR, Failed to create molecule", input_count, enumerated_count, smi)
                        continue

                    molh = Chem.AddHs(mol)

                    mol_with_confs = gen_conformers(molh, rms_threshold, minimize_cycles, True, num_conformers=num_conformers)
                    for idx in range(mol_with_confs.GetNumConformers()):
                        # mol_with_confs.SetDoubleProp('Energy',
                        #                              mol_with_confs.GetConformer(idx).GetDoubleProp('Energy'))
                        # mol_with_confs.SetDoubleProp('Energy_Delta',
                        #                              mol_with_confs.GetConformer(idx).GetDoubleProp(
                        #                                  'Energy_Delta'))
                        cxsmi = Chem.MolToCXSmiles(mol_with_confs)
                        writer.write("{}\t{}\t{}\t{}\t{}\n".format(cxsmi, id, str(idx),
                                                         mol_with_confs.GetConformer(idx).GetDoubleProp('Energy'),
                                                         mol_with_confs.GetConformer(idx).GetDoubleProp( 'Energy_Delta')))
                        conformer_count += 1
                        conf_count_for_mol += 1

                except KeyboardInterrupt:
                    utils.log('Interrupted')
                    return input_count, enumerated_count, conformer_count, error_count
                except:
                    error_count += 1
                    traceback.print_exc()

                utils.log('INFO, Generated', conf_count_for_mol, 'conformers for', smi)

    return input_count, enumerated_count, conformer_count, error_count


### start main execution #########################################

def main():

    # Example:
    #   python -m moldb.conformers.py -i foo.smi -o bar.cxsmi

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate conformers')
    parser.add_argument('-i', '--input', required=True, help="Input file as SMILES")
    parser.add_argument('-o', '--output', required=True, help="Output file as CXSMILES")
    parser.add_argument('-n', '--num-conformers', type=int,
                        help="Number of conformers to generate. If not specified the Inhibox rules are used")
    parser.add_argument('-m', '--minimize-cycles', type=int, default=500, help="Number of MMFF minimisation cycles")
    parser.add_argument('-t', '--rms-threshold', type=float, default=1.0, help="RMS threshold for excluding conformers")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("le_conformers: ", args)

    start = time.time()
    input_count, enumerated_count, conformer_count, error_count = \
        execute(args.input, args.output,
                num_conformers=args.num_conformers,
                minimize_cycles=args.minimize_cycles,
                rms_threshold=args.rms_threshold,
                interval=args.interval
                )
    end = time.time()

    DmLog.emit_event('Inputs:', input_count, 'Enumerated:', enumerated_count,
                     'Conformers:', conformer_count, 'Errors:', error_count, 'Time (s):', end - start)
    DmLog.emit_cost(conformer_count)


if __name__ == "__main__":
    main()
