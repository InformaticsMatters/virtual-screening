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

from oddt import toolkit, fingerprints

fingerprint_types = {
    'if': fingerprints.InteractionFingerprint,
    'sif': fingerprints.SimpleInteractionFingerprint,
    'splif': fingerprints.SPLIF
}

metric_types = {
    'tanimoto': fingerprints.tanimoto,
    'dice': fingerprints.dice,
    'splif': fingerprints.similarity_SPLIF
}


def execute(inputs, protein, ligands, output, fingerprint='if', metric='dice', field_name='ODDT_IFP', interval=0):

    count = 0
    success = 0
    errors = 0

    utils.expand_path(output)

    if fingerprint == 'splif':
        metric = 'splif'

    protein_mol = next(toolkit.readfile('pdb', protein))
    protein_mol.protein = True
    ligand_mols = toolkit.readfile('sdf', ligands)
    input_mols = toolkit.readfile('sdf', inputs)
    fps = [ fingerprint_types[fingerprint](m, protein_mol) for m in ligand_mols ]

    writer = toolkit.Outputfile('sdf', output, overwrite=True)

    try:
        for mol in input_mols:
            count += 1
            if not mol:
                DmLog.emit_event("Failed to read molecule", count)
                errors += 1
                continue
            try:
                success += 1
                fp2 = fingerprint_types[fingerprint](mol, protein_mol)
                best = 0
                sum = 0
                for fp in fps:
                    score = metric_types[metric](fp, fp2)
                    sum += score
                    if score > best:
                        best = score
                mol.data[field_name + '_MAX'] = best
                mol.data[field_name + '_AVG'] = sum / len(fps)
                writer.write(mol)

            except RuntimeError as e:
                utils.log("Failed", count, e)
                errors += 1

    finally:
        writer.close()

        if interval and count % interval == 0:
            DmLog.emit_event("Processed {} records, {} errors".format(count, errors))
        if success % 10000 == 0:
            DmLog.emit_cost(success)

    return count, errors


def main():

    # Example:
    #   ./oddt_ifp.py -i data/dhfr_candidates.sdf -p data/dhfr-receptor.pdb -l data/dhfr-ligand.mol -o foo.sdf

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Minimize structures')
    parser.add_argument('-i', '--input', required=True, help="File with ligands to score (.sdf)")
    parser.add_argument('-p', '--protein', required=True, help="Protein in PDB format")
    parser.add_argument('-l', '--ligands', required=True, help="Known ligands (.sdf)")
    parser.add_argument('-o', '--output', required=True, help="Output file (.sdf)")
    parser.add_argument('-f', '--fingerprint', default='if', choices=['if', 'sif', 'splif'],
                        help="Which fingerprint (if=InteractionFingerprint), sif=SimpleInteractionFingerprint" +
                        " splif=SPLIF")
    parser.add_argument('-m', '--metric', default='dice', choices=['dice', 'tanimoto'],
                        help="Which metric (ignored if using splif fingerprint")
    parser.add_argument('--field-name', default='ODDT_IFP', help="Base name for the output fields")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log("minimize.py: ", args)

    count, errors = execute(args.input, args.protein, args.ligands, args.output,
                            fingerprint=args.fingerprint, metric=args.metric,
                            field_name=args.field_name, interval=args.interval)
    DmLog.emit_event('Processed {} molecules. {} errors'.format(count, errors))
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()

