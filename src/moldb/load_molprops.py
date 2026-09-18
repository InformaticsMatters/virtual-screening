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


def load_molprops(input_files):

    with Session(engine) as session:
        cursor = session.connection().connection.cursor()

        for input_file in input_files:
            letters = string.ascii_lowercase
            table_name = 'props_' + ''.join(random.choice(letters) for i in range(16))
            sql = "COPY {}(smiles, id, hac, rot_bonds, ring_count, aromatic_ring_count, chiral_centres, undef_chiral_centres, num_sp3, logp, tpsa) FROM STDIN DELIMITER E'\t' CSV;".format(table_name)
            DmLog.emit_event('Copying properties using:\n', sql)

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

            q = "UPDATE molecule m SET hac=c.hac, rot_bonds=c.rot_bonds, ring_count=c.ring_count, \n" + \
                "  aromatic_ring_count=c.aromatic_ring_count, chiral_centres=c.chiral_centres, \n" + \
                "  undef_chiral_centres=c.undef_chiral_centres, num_sp3=c.num_sp3, logp=c.logp, tpsa=c.tpsa\n" + \
                "  FROM {} as c WHERE m.id = c.id"
            stmt = text(q.format(table_name))
            session.connection().execute(stmt)

            stmt = text("DROP TABLE {}".format(table_name))
            session.connection().execute(stmt)

            session.commit()


def main():

    ### typical usage:
    # python -m moldb.load_molprops --inputs inputs.smi

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Load molecular properties')
    parser.add_argument("--inputs", nargs='+', required=True, help="Input files to load")

    args = parser.parse_args()
    DmLog.emit_event("load_molprops: ", args)

    t0 = time.time()
    load_molprops(args.inputs)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processing complete in {} seconds'.format(duration_s))
    # DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
