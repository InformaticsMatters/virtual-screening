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

import os, argparse, time, random, string, io

from sqlalchemy import MetaData, Table, Column, Text
from sqlalchemy import text
from sqlalchemy.orm import sessionmaker, Session

from . import models
from . models import Library, File

from . import standardize
from dm_job_utilities.dm_log import DmLog

engine = models.get_engine(echo=False)


def run1(input, library_name, delimiter, name_column=None, skip_lines=None, interval=None):

    # 1. standardize the molecules
    standardized_file, count, errors = _standardize(input, delimiter,
                                                    name_column=name_column, skip_lines=skip_lines, interval=interval)

    # 2. load into temp db table
    table_name = _load_standardized(standardized_file)

    # 3. Insert library, filename and do the bulk transfer into molecules
    _load_molecules(library_name, input, table_name)

    return count, errors


def _standardize(input, delimiter, name_column=None, skip_lines=None, interval=None):

    DmLog.emit_event('Standardizing', input)

    input_file = os.path.basename(input)
    base = ['outputs', 'standardize']
    if not os.path.isdir(os.path.join(*base)):
        os.makedirs(os.path.join(*base), exist_ok=True)
    path = base
    path.append(input_file)
    output_file = os.path.join(*path)

    count, errors = standardize.standardize([input], output_file, delimiter,
                                            name_column=name_column, skip_lines=skip_lines, interval=interval)
    return output_file, count, errors


def _load_standardized(molecules):

    DmLog.emit_event('Loading standardized molecules:', molecules)

    letters = string.ascii_lowercase
    table_name = 'std_' + ''.join(random.choice(letters) for i in range(16))
    sql = "COPY {}(smiles, orig_smiles, id) FROM STDIN DELIMITER E'\t' CSV HEADER;".format(table_name)

    session = sessionmaker(bind=engine)()
    cursor = session.connection().connection.cursor()

    metadata_obj = MetaData()
    table = Table(table_name, metadata_obj,
                  Column('smiles', Text),
                  Column('orig_smiles', Text),
                  Column('id', Text)
                  )
    metadata_obj.create_all(engine)

    input = open(molecules, 'r')
    cursor.copy_expert(sql, input)

    session.commit()

    return table_name


def _load_molecules(library_name, file_name, table_name):

    DmLog.emit_event('Loading', file_name, 'into', table_name)

    with Session(engine) as session:

        query = session.query(Library).filter_by(name=library_name)
        l = query.first()
        if not l:
            l = Library(name=library_name)
            session.add(l)

        f = File(name=file_name, library_id=l.id, load_table=table_name)
        l.files.append(f)

        session.flush()

        # pull the smiles into the molecule table
        q = "INSERT INTO molecule(smiles, created_at) SELECT smiles, now() FROM {}" +\
            " ON CONFLICT ON CONSTRAINT molecule_smiles_key DO NOTHING"
        stmt = text(q.format(table_name))
        session.connection().execute(stmt)

        # pull the data into the supply table
        q = "WITH rec AS (" +\
            "  SELECT i.orig_smiles, i.id AS code, m.id AS molid FROM {} i JOIN molecule m on m.smiles = i.smiles)\n" +\
            "  INSERT INTO supply(created_at, smiles, code, molecule_id, file_id)\n" +\
            "  SELECT now(), orig_smiles, code, molid, :f FROM rec"
        stmt = text(q.format(table_name))
        session.connection().execute(stmt, {'f': f.id})

        session.commit()

        return l.id, f.id


def main():

    ### typical usage:
    # python -m moldb.run1-load-library -i data/10000.smi -l chemspace -n 1 --interval 1000

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Fist run stage of compound prep')
    parser.add_argument('-i', '--input', help="Input file as SMILES")
    parser.add_argument('-l', '--library', required=True, help="Library name")
    parser.add_argument('-d', '--delimiter', default='\t', help="Delimiter")
    parser.add_argument('-n', '--name-column', required=True, type=int, help="Column for name field (zero based)")
    parser.add_argument('--skip-lines', default=0, type=int, help="Skip this many lines e.g. use 1 to skip header line")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("run1-load-library: ", args)

    t0 = time.time()
    count, errors = run1(args.input, args.library, args.delimiter, name_column=args.name_column, interval=args.interval,
                         skip_lines=args.skip_lines)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processed {} records in {} seconds, {} errors.'.format(count, duration_s, errors))
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()

