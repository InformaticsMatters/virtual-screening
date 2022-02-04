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

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign


def minimise(input, output, cycles=200, interval=0):

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
                    utils.log_dm_event("Failed to read molecule", count)
                    errors += 1
                    continue
                try:
                    mol = Chem.RemoveHs(mol)
                    molh = Chem.AddHs(mol, addCoords=True)
                    res = AllChem.MMFFOptimizeMoleculeConfs(molh, maxIters=cycles)
                    converged, energy = res[0]
                    if converged == 1:
                        non_converged += 1
                    elif converged == -1:
                        utils.log_dm_event("Force field could not be set up for molecule", count)
                        errors += 1

                    probe_mol = Chem.RemoveHs(molh)

                    rmsd = rdMolAlign.AlignMol(probe_mol, mol)
                    probe_mol.SetDoubleProp('RMSD', rmsd)
                    probe_mol.SetIntProp('CONVERGED', converged)
                    probe_mol.SetDoubleProp('ENERGY', energy)

                    writer.write(probe_mol)
                    success += 1

                except RuntimeError as e:
                    utils.log("Failed to minimize", count, Chem.MolToSmiles(mol), e)
                    errors += 1

                if interval and count % interval == 0:
                    utils.log_dm_event("Processed {} records, {} errors".format(count, errors))
                if success % 10000 == 0:
                    utils.log_dm_cost(success)

    return count, success, errors, non_converged


def main():

    # Example:
    #   ./minimize.py -i input.sdf -o output.sdf

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Minimize structures')
    parser.add_argument('-i', '--input', required=True, help="File with inputs")
    parser.add_argument('-o', '--output', required=True, help="Output file")
    parser.add_argument('-c', '--cycles', type=int, default=200, help="Number of minimization cycles")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log("minimize.py: ", args)

    count, success, errors, non_converged = minimise(args.input, args.output, args.cycles, interval=args.interval)
    utils.log_dm_event('Processed {} molecules. {} did not converge. {} errors'.format(count, non_converged, errors))
    utils.log_dm_cost(success)
    
    
if __name__ == "__main__":
    main()
        
