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

import argparse, time, os, random, string, traceback

from . import models

from sqlalchemy import text
from sqlalchemy.orm import Session
from sqlalchemy import MetaData, Table, Column, Integer, SmallInteger, String, Text, Float

from dm_job_utilities.dm_log import DmLog

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen

engine = models.get_engine(echo=False)


def calc_molprops(input_file, output_file, interval=None):

    DmLog.emit_event('Writing calculated properties to', output_file)

    count = 0
    errors = 0

    with open(output_file, 'wt') as writer:
        with open(input_file, 'rt') as reader:
            for line in reader:
                count += 1
                if interval and count % interval == 0:
                    DmLog.emit_event("Processed {} records".format(count))

                tokens = line.strip().split('\t')
                mol = Chem.MolFromSmiles(tokens[0])
                values = []
                values.append(tokens[0])
                values.append(tokens[1])
                if not mol:
                    errors += 1
                    DmLog.emit_event("Failed to process record", count)
                    continue

                try:
                    values.append(str(mol.GetNumHeavyAtoms()))
                    values.append(str(rdMolDescriptors.CalcNumRotatableBonds(mol)))
                    values.append(str(rdMolDescriptors.CalcNumRings(mol)))
                    values.append(str(rdMolDescriptors.CalcNumAromaticRings(mol)))
                    num_cc, num_undef_cc = _get_num_chiral_centers(mol)
                    values.append(str(num_cc))
                    values.append(str(num_undef_cc))
                    values.append(str(sum((x.GetHybridization() == Chem.HybridizationType.SP3) for x in mol.GetAtoms())))
                    logp = Crippen.MolLogP(mol)
                    # logp values have a silly number of decimal places
                    values.append("%.2f" % logp)
                    tpsa = rdMolDescriptors.CalcTPSA(mol)
                    # tpsa values have a silly number of decimal places
                    values.append("%.2f" % tpsa)

                    writer.write("\t".join(values) + '\n')

                except:
                   errors += 1
                   DmLog.emit_event('Failed to process record', count)
                   #traceback.print_exc()
                   continue

    DmLog.emit_event('Calculated properties for', count, 'molecules.', errors, 'errors')

    return output_file, count, errors


def _get_num_chiral_centers(mol):
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    undef_cc = 0
    for cc in chiral_centers:
        if cc[1] == '?':
            undef_cc += 1
    return len(chiral_centers), undef_cc


def main():

    ### typical usage:
    # python -m moldb.calc_molprops --count 1000 --interval 100

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Calculate molecular properties')
    parser.add_argument('-i', '--input', required=True, help="Input file (.smi)")
    parser.add_argument('-o', '--output', required=True, help="Output file (.smi)")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("calc_molprops: ", args)

    t0 = time.time()
    calc_molprops(args.input, args.output, interval=args.interval)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processing complete in {} seconds'.format(duration_s))
    # DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
