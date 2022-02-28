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


import sys, argparse, time
import utils
from dm_job_utilities.dm_log import DmLog

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdForceFieldHelpers


def minimize(input, output, cycles=200, molecule='minimized', interval=0):

    count = 0
    errors = 0
    success = 0
    non_converged = 0

    utils.expand_path(output)

    with Chem.SDMolSupplier(input) as suppl:
        with Chem.SDWriter(output) as writer:
            for mol in suppl:
                count += 1
                if not mol:
                    DmLog.emit_event("Failed to read molecule", count)
                    errors += 1
                    continue
                try:
                    mol = Chem.RemoveHs(mol)
                    molh = Chem.AddHs(mol, addCoords=True)

                    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(molh, AllChem.MMFFGetMoleculeProperties(molh))
                    energy0 = ff.CalcEnergy()

                    res = AllChem.MMFFOptimizeMoleculeConfs(molh, maxIters=cycles)
                    converged, energy1 = res[0]
                    if converged == 1:
                        non_converged += 1
                    elif converged == -1:
                        DmLog.emit_event("Force field could not be set up for molecule", count)
                        errors += 1

                    probe_mol = Chem.RemoveHs(molh)
                    rmsd = rdMolAlign.AlignMol(probe_mol, mol)

                    if molecule == 'original':
                        w_mol = mol
                    elif molecule == 'minimized':
                        w_mol = probe_mol
                    elif molecule == 'merged':
                        w_mol = Chem.RWMol(mol)
                        w_mol.InsertMol(probe_mol)
                        Chem.SanitizeMol(w_mol)
                    else:
                        raise ValueError('Invalid molecule to write:', molecule)

                    w_mol.SetDoubleProp('RMSD', rmsd)
                    if converged:
                        converged_b = 'no'
                    else:
                        converged_b = 'yes'
                    w_mol.SetProp('CONVERGED', converged_b)
                    w_mol.SetDoubleProp('ENERGY_MIN', energy1)
                    w_mol.SetDoubleProp('ENERGY_DELTA', energy0 - energy1)

                    writer.write(w_mol)
                    success += 1

                except RuntimeError as e:
                    utils.log("Failed to minimize", count, Chem.MolToSmiles(mol), e)
                    errors += 1

                if interval and count % interval == 0:
                    DmLog.emit_event("Processed {} records, {} errors".format(count, errors))
                if success % 10000 == 0:
                    DmLog.emit_cost(success)

    return count, success, errors, non_converged


def main():

    # Example:
    #   ./minimize.py -i input.sdf -o output.sdf

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Minimize structures')
    parser.add_argument('-i', '--input', required=True, help="File with inputs")
    parser.add_argument('-o', '--output', required=True, help="Output file")
    parser.add_argument('-c', '--cycles', type=int, default=200, help="Number of minimization cycles")
    parser.add_argument('-m', '--molecule', choices=['minimized', 'original', 'merged'], default='minimized',
                        help='Which molecule to write in the output SD file')
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log("minimize.py: ", args)

    count, success, errors, non_converged = minimize(args.input, args.output, args.cycles,
                                                     molecule=args.molecule, interval=args.interval)
    DmLog.emit_event('Processed {} molecules. {} did not converge. {} errors'.format(count, non_converged, errors))
    DmLog.emit_cost(success)
    
    
if __name__ == "__main__":
    main()
        
