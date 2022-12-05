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

import argparse, datetime, time

import utils
from . import models, filter, moldb_utils

from dm_job_utilities.dm_log import DmLog

from sqlalchemy import text



engine = models.get_engine(echo=False)


def analyse(outfile, specification, skip_enum, skip_conf):
    utils.expand_path(outfile)
    DmLog.emit_event('Writing to', outfile)

    if specification:
        filters = moldb_utils.read_specification(specification)
    else:
        filters = {}

    with open(outfile, 'wt') as w:
        analyseEnvironment(w)
        if filters:
            w.write("Specification: {}\n".format(specification))
            f_sql = filter._gen_filters(filters)
            w.write('Filter terms: {}\nFilter SQL: {}\n'.format(str(filters), f_sql))
        analyseTables(w, filters, skip_enum, skip_conf)


def analyseEnvironment(writer):

    url = models.gen_url(obscure_password=True)
    writer.write('Database: {}\n'.format(url))
    writer.write('Date: {}\n'.format(datetime.datetime.now()))


def analyseTables(writer, filters, skip_enum, skip_conf):
    writer.write('Database analysis:\n')
    with engine.connect() as conn:
        analyseLibrary(writer, filters, conn)
        analyseMolecule(writer, filters, conn)
        if not skip_enum:
            analyseEnumeration(writer, filters, conn)
        if not skip_conf:
            analyseConformer(writer, filters, conn)


def analyseLibrary(writer, filters, conn):
    DmLog.emit_event('Analysing library')
    sql = 'SELECT name FROM library'
    t0 = time.time()
    results = conn.execute(text(sql))
    names = []
    for row in results:
        lib_name = row[0]
        names.append(lib_name)
    t1 = time.time()
    msg = '  library: {}\n'.format(', '.join(names))
    msg += '           took {}s\n'.format((t1-t0))

    writer.write(msg)
    writer.flush()


def analyseMolecule(writer, filters, conn):
    DmLog.emit_event('Analysing molecule')

    t0 = time.time()

    sql = 'SELECT COUNT(*) FROM molecule'
    results = conn.execute(text(sql))
    row = results.fetchone()
    mol_count = row[0]

    sql = 'SELECT COUNT(*) FROM molecule WHERE hac IS NULL'
    results = conn.execute(text(sql))
    row = results.fetchone()
    need_calc = row[0]

    if filters:
        f_sql = filter._gen_filters(filters)
        DmLog.emit_event('Using filter:', f_sql)
        sql = 'SELECT COUNT(*) FROM molecule WHERE ' + f_sql
        results = conn.execute(text(sql))
        row = results.fetchone()
        filt_count = row[0]
        msg = '  molecule: {} rows, {} pass filters, {} need property calculation\n'.format(mol_count, filt_count, need_calc)
    else:
        msg = '  molecule: {} rows, {} need property calculation\n'.format(mol_count, need_calc)

    t1 = time.time()
    msg += '            took {}s\n'.format((t1-t0))

    writer.write(msg)
    writer.flush()


def analyseEnumeration(writer, filters, conn):
    DmLog.emit_event('Analysing enumeration')

    t0 = time.time()

    sql = 'SELECT COUNT(*) FROM enumeration'
    results = conn.execute(text(sql))
    row = results.fetchone()
    enum_count = row[0]

    sql = 'SELECT count(DISTINCT molecule_id) FROM enumeration'
    results = conn.execute(text(sql))
    row = results.fetchone()
    enum_mols = row[0]

    sql = 'SELECT COUNT(*) FROM molecule WHERE NOT EXISTS (SELECT 1 FROM enumeration WHERE enumeration.molecule_id = molecule.id)'
    results = conn.execute(text(sql))
    row = results.fetchone()
    need_enum = row[0]

    msg = '  enumeration - all data: {} rows, {} molecules enumerated, {} need enumeration\n'.format(enum_count, enum_mols, need_enum)

    if filters:
        f_sql = filter._gen_filters(filters, prefix='m.')

        sql = 'SELECT COUNT(*) FROM enumeration e JOIN molecule m ON m.id = e.molecule_id WHERE' + f_sql
        results = conn.execute(text(sql))
        row = results.fetchone()
        enum_filt = row[0]

        sql = 'SELECT count(DISTINCT molecule_id) FROM enumeration e JOIN molecule m ON m.id = e.molecule_id WHERE' + f_sql
        results = conn.execute(text(sql))
        row = results.fetchone()
        mols_filt = row[0]

        sql = 'SELECT COUNT(*) FROM molecule m WHERE NOT EXISTS (SELECT 1 FROM enumeration e WHERE e.molecule_id = m.id) AND' + f_sql
        results = conn.execute(text(sql))
        row = results.fetchone()
        need_enum_filt = row[0]

        msg += '             with filter: {} rows, {} molecules enumerated, {} need enumeration\n'.format(enum_filt, mols_filt, need_enum_filt)

    t1 = time.time()

    msg += '             took {}s\n'.format((t1-t0))

    writer.write(msg)
    writer.flush()


def analyseConformer(writer, filters, conn):
    DmLog.emit_event('Analysing conformer')

    t0 = time.time()

    sql = 'SELECT count(*) FROM conformer'
    results = conn.execute(text(sql))
    row = results.fetchone()
    conf_count = row[0]

    # sql = 'SELECT count(DISTINCT molecule_id) FROM enumeration WHERE EXISTS (SELECT 1 FROM conformer WHERE conformer.enumeration_id = enumeration.id)'
    # results = conn.execute(text(sql))
    # row = results.fetchone()
    # enum_mols = row[0]

    # sql = 'SELECT count(DISTINCT molecule_id) FROM enumeration WHERE NOT EXISTS (SELECT 1 FROM conformer WHERE conformer.enumeration_id = enumeration.id)'
    # results = conn.execute(text(sql))
    # row = results.fetchone()
    # mols_need_conf = row[0]

    sql = 'SELECT count(*) FROM enumeration WHERE NOT EXISTS (SELECT 1 FROM conformer WHERE conformer.enumeration_id = enumeration.id)'
    results = conn.execute(text(sql))
    row = results.fetchone()
    enums_need_conf = row[0]

    # msg = '  conformer - all data: {} rows, {} molecules have conformers, {} molecules need conformer generation, {} enumerated forms need conformer generation\n'\
    #     .format(conf_count, enum_mols, mols_need_conf, enums_need_conf)
    msg = '  conformer - all data: {} rows, {} enumerations need conformer generation\n' \
        .format(conf_count, enums_need_conf)

    if filters:
        f_sql = filter._gen_filters(filters, prefix='m.')

        sql = 'SELECT count(*) FROM conformer c JOIN enumeration e ON e.id = c.enumeration_id JOIN molecule m ON m.id = e.molecule_id WHERE' + f_sql
        results = conn.execute(text(sql))
        row = results.fetchone()
        filt_count = row[0]

        # sql = 'SELECT count(DISTINCT molecule_id) FROM enumeration e JOIN molecule m ON e.molecule_id = m.id ' + \
        #       'WHERE EXISTS (SELECT 1 FROM conformer WHERE conformer.enumeration_id = e.id) AND' + f_sql
        # results = conn.execute(text(sql))
        # row = results.fetchone()
        # enum_mols_filt = row[0]

        # sql = 'SELECT count(DISTINCT molecule_id) FROM enumeration e JOIN molecule m ON e.molecule_id = m.id' + \
        #       ' WHERE NOT EXISTS (SELECT 1 FROM conformer WHERE conformer.enumeration_id = e.id) AND' + f_sql
        # results = conn.execute(text(sql))
        # row = results.fetchone()
        # mols_need_conf_filt = row[0]

        sql = 'SELECT count(*) FROM enumeration e JOIN molecule m ON e.molecule_id = m.id ' + \
              'WHERE NOT EXISTS (SELECT 1 FROM conformer WHERE conformer.enumeration_id = e.id) AND' + f_sql
        results = conn.execute(text(sql))
        row = results.fetchone()
        enums_need_conf_filt = row[0]

        # msg += '             with filter: {} rows, {} molecules have conformers, {} molecules need conformer generation, {} enumerated forms need conformer generation\n'\
        #     .format(filt_count, enum_mols_filt, mols_need_conf_filt, enums_need_conf_filt)
        msg += '           with filter: {} rows, {} enumerations need conformer generation\n' \
            .format(filt_count, enums_need_conf_filt)

    t1 = time.time()
    msg += '           took {}s\n'.format((t1-t0))

    writer.write(msg)


def main():

    # Example:
    #   python -m moldb.analyse -s specification.txt -o report.txt

    # ----- command line args definitions ---------------------------------------------

    parser = argparse.ArgumentParser(description='Analyse')
    parser.add_argument('-s', '--specification', help="Filter specification file")
    parser.add_argument('-o', '--output', default='report.txt', help="Output file")
    parser.add_argument('--skip-enumeration', action='store_true', help="Don't analyse the enumeration table")
    parser.add_argument('--skip-conformer', action='store_true', help="Don't analyse the conformer table")

    args = parser.parse_args()
    DmLog.emit_event("analyse: ", args)

    if args.specification:
        filters = moldb_utils.read_specification(args.specification)
    else:
        filters = {}

    analyse(args.output, args.specification, args.skip_enumeration, args.skip_conformer)


if __name__ == "__main__":
    main()
