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
Loads standardize molecules.
"""

import argparse, string, random


from sqlalchemy import MetaData, Table, Column, Text
from sqlalchemy import text
from sqlalchemy.orm import sessionmaker, Session

from . import models
from . models import Library, File

from dm_job_utilities.dm_log import DmLog


engine = models.get_engine(echo=False)


def load_data(inputs, library_name):
    """
    Load standardized molecules into the database.

    :param inputs: file name to load (contents: std_smi\torig_smi\tid
    :param library_name: Name of the library
    :return:
    """

    session = sessionmaker(bind=engine)()
    cursor = session.connection().connection.cursor()

    query = session.query(Library).filter_by(name=library_name)
    lib = query.first()
    if not lib:
        lib = Library(name=library_name)
        session.add(lib)
        session.flush()

    for input in inputs:
        DmLog.emit_event('Loading standardized molecules:', input)

        letters = string.ascii_lowercase
        table_name = 'std_' + ''.join(random.choice(letters) for i in range(16))
        sql = "COPY {}(smiles, orig_smiles, id) FROM STDIN DELIMITER E'\t' CSV HEADER;".format(table_name)

        metadata_obj = MetaData()
        tmp_table = Table(table_name, metadata_obj,
                      Column('smiles', Text),
                      Column('orig_smiles', Text),
                      Column('id', Text)
                      )
        metadata_obj.create_all(engine)

        reader = open(input, 'r')
        cursor.copy_expert(sql, reader)

        # now transfer the data into our tables
        file = File(name=input, library_id=lib.id, load_table=table_name)
        lib.files.append(file)

        session.flush()

        # pull the smiles into the molecule table
        q = "INSERT INTO molecule(smiles, created_at) SELECT smiles, now() FROM {}" + \
            " ON CONFLICT ON CONSTRAINT molecule_smiles_key DO NOTHING"
        stmt = text(q.format(table_name))
        session.connection().execute(stmt)

        # pull the data into the supply table
        q = "WITH rec AS (" + \
            "  SELECT i.orig_smiles, i.id AS code, m.id AS molid FROM {} i JOIN molecule m on m.smiles = i.smiles)\n" + \
            "  INSERT INTO supply(created_at, smiles, code, molecule_id, file_id)\n" + \
            "  SELECT now(), orig_smiles, code, molid, :f FROM rec"
        stmt = text(q.format(table_name))
        session.connection().execute(stmt, {'f': file.id})

        stmt = text("DROP TABLE {}".format(table_name))
        session.connection().execute(stmt)

        session.commit()


def main():

    # Example:
    #   python -m moldb.load_standardized -i standardized.smi

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='load_enums')
    parser.add_argument('-i', '--inputs', nargs='+', help='Input file (.smi)')
    parser.add_argument('-l', '--library-name', required=True, help='Library name')

    args = parser.parse_args()
    DmLog.emit_event("load_standardized: ", args)

    load_data(args.inputs, args.library_name)


if __name__ == "__main__":
    main()
