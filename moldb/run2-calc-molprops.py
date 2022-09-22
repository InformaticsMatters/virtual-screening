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


def run2(count=None, interval=None):

    # 1. extract molecules needing calculation
    file_name = _extract_mols(count=count)

    # 2. calculate properties
    output_file, count, errors = _calc_molprops(file_name, interval=interval)

    # 3. load properties
    _load_molprops(output_file)


def _extract_mols(count=None):

    path = ['outputs', 'calculate']
    path_s = os.path.join(*path)
    if not os.path.isdir(path_s):
        os.makedirs(path_s, exist_ok=True)

    path.append('mols_in_' + ''.join(random.choice(string.ascii_lowercase) for i in range(16)))
    file_name = os.path.join(*path)

    DmLog.emit_event('Writing to', file_name)

    with Session(engine) as session:

        if count:
            sql = "COPY (SELECT smiles, id FROM molecule WHERE hac IS NULL LIMIT {}) TO STDOUT".format(count)
        else:
            sql = "COPY (SELECT smiles, id FROM molecule WHERE hac IS NULL) TO STDOUT"

        DmLog.emit_event('Extracting uncalculated molecules using:\n', sql)

        output = open(file_name, 'w')
        cursor = session.connection().connection.cursor()
        cursor.copy_expert(sql, output)

        session.commit()

    return file_name


def _calc_molprops(file_name, interval=None):


    output_file = file_name.replace('/mols_in_', '/mols_out_', 1)
    DmLog.emit_event('Writing calculated properties to', output_file)

    count = 0
    errors = 0

    with open(output_file, 'wt') as writer:
        with open(file_name, 'rt') as reader:
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


def _load_molprops(input_file):

    letters = string.ascii_lowercase
    table_name = 'props_' + ''.join(random.choice(letters) for i in range(16))
    sql = "COPY {}(smiles, id, hac, rot_bonds, ring_count, aromatic_ring_count, chiral_centres, undef_chiral_centres, num_sp3, logp, tpsa) FROM STDIN DELIMITER E'\t' CSV HEADER;".format(table_name)
    DmLog.emit_event('Copying properties using:\n', sql)

    with Session(engine) as session:
        cursor = session.connection().connection.cursor()

        metadata_obj = MetaData()
        table = Table(table_name, metadata_obj,
                    Column('id', Integer),
                    Column('smiles', Text),
                    Column('hac', SmallInteger),
                    Column('rot_bonds', SmallInteger),
                    Column('ring_count', SmallInteger),
                    Column('aromatic_ring_count', SmallInteger),
                    Column('chiral_centres', SmallInteger),
                    Column('undef_chiral_centres', SmallInteger),
                    Column('num_sp3', SmallInteger()),
                    Column('logp', Float),
                    Column('tpsa', Float)
                      )
        metadata_obj.create_all(engine)

        input = open(input_file, 'r')
        cursor.copy_expert(sql, input)

        q = "UPDATE molecule m SET hac=c.hac, rot_bonds=c.rot_bonds, ring_count=c.ring_count, \n" +\
            "  aromatic_ring_count=c.aromatic_ring_count, chiral_centres=c.chiral_centres, \n" +\
            "  undef_chiral_centres=c.undef_chiral_centres, num_sp3=c.num_sp3, logp=c.logp, tpsa=c.tpsa\n" +\
            "  FROM {} as c WHERE m.id = c.id"
        stmt = text(q.format(table_name))
        session.connection().execute(stmt)

        session.commit()


def _get_num_chiral_centers(mol):
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    undef_cc = 0
    for cc in chiral_centers:
        if cc[1] == '?':
            undef_cc += 1
    return len(chiral_centers), undef_cc


def main():

    ### typical usage:
    # python -m moldb.run2-calc-molprops --count 1000 --interval 100

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Second run stage of compound prep')
    parser.add_argument("--count", type=int, help="Max number of molecules to work on")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("run2-calc-molprops: ", args)

    t0 = time.time()
    run2(count=args.count, interval=args.interval)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processing complete in {} seconds'.format(duration_s))
    # DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
