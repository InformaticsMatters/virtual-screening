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


import os, glob, argparse, time
import utils

from sqlalchemy import text
from sqlalchemy.orm import Session

from . import models, moldb_utils

import rdkit_utils

from rdkit import Chem

from dm_job_utilities.dm_log import DmLog



engine = models.get_engine(echo=False)


def _add_term(sql, min_value, max_value, column, count, prefix=''):

    if min_value is None and max_value is None:
        return sql, count

    if min_value is not None and max_value is not None and min_value == max_value:
        if count:
            sql = sql + ' AND '
        sql = sql + prefix + column + ' = ' + str(min_value)
        count += 1
    else:
        if min_value is not None:
            float(min_value)
            if count:
                sql = sql + ' AND '
            sql = sql + prefix + column + ' >= ' + str(min_value)
            count += 1
        if max_value is not None:
            float(max_value)
            if count:
                sql = sql + ' AND '
            sql = sql + prefix + column + ' <= ' + str(max_value)
            count += 1
    return sql, count


def _gen_filters(filters, smarts, prefix=''):

    terms = (
        (filters.get('min_hac'), filters.get('max_hac'), 'hac'),
        (filters.get('min_rotb'), filters.get('max_rotb'), 'rot_bonds'),
        (filters.get('min_rings'), filters.get('max_rings'), 'ring_count'),
        (filters.get('min_aro_rings'), filters.get('max_aro_rings'), 'aromatic_ring_count'),
        (filters.get('min_chiral_centres'), filters.get('max_chiral_centres'), 'chiral_centres'),
        (filters.get('min_undefined_chiral_centres'), filters.get('max_undefined_chiral_centres'), 'undef_chiral_centres'),
        (filters.get('min_sp3'), filters.get('max_sp3'), 'num_sp3'),
        (filters.get('min_logp'), filters.get('max_logp'), 'logp'),
        (filters.get('min_tpsa'), filters.get('max_tpsa'), 'tpsa')
    )
    sql = " "
    count = 0
    for term in terms:
        sql, count = _add_term(sql, term[0], term[1], term[2], count, prefix=prefix)

    if smarts:
        if count:
            sql += ' AND ('
            sql = add_smarts_queries(sql, smarts, prefix)
        if count:
            sql += ')'

    return sql


def add_smarts_queries(sql, smarts, prefix):
    if smarts:
        for i, s in enumerate(smarts):
            if i:
                sql += ' OR '
            sql += prefix + "mol @> '" + s + "'::qmol"
    return sql


def _do_filter(sql, output_file, filters, smarts, dry_run=False, prefix='', suffix=''):

    utils.expand_path(output_file)

    filter_sql = _gen_filters(filters, smarts, prefix=prefix)

    sql = sql + filter_sql + suffix

    utils.log('Query SQL:\n', sql)

    if not dry_run:
        with Session(engine) as session:

            DmLog.emit_event('Executing query')

            output = open(output_file, 'w')
            cursor = session.connection().connection.cursor()
            cursor.copy_expert(sql, output)

            session.commit()


def filter_molecules(output_file, count, filters, smarts, dry_run=False):

    sql = "COPY (SELECT m.smiles, m.id, " +\
          "concat_ws(',', VARIADIC array_agg(s.code)) codes, " +\
          "concat_ws(',', VARIADIC array_agg(l.name)) libs  " +\
          "FROM molecule m " +\
          "JOIN supply s on s.molecule_id = m.id " +\
          "JOIN file f ON f.id = s.file_id " +\
          "JOIN library l ON l.id = f.library_id " +\
          "GROUP BY m.id HAVING"
    if count:
        suffix = ' LIMIT ' + str(int(count)) + ') TO STDOUT'
    else:
        suffix = ') TO STDOUT'

    _do_filter(sql, output_file, filters, smarts, dry_run=dry_run, prefix='m.', suffix=suffix)


def filter_need_enum(output_file, count, filters, smarts, dry_run=False):

    utils.expand_path(output_file)

    sql = "COPY (SELECT m.smiles, m.id FROM molecule m WHERE"
    if count:
        suffix = ' AND NOT EXISTS (SELECT 1 FROM enumeration e WHERE e.molecule_id = m.id) LIMIT ' + \
                 str(int(count)) + ') TO STDOUT'
    else:
        suffix = ' AND NOT EXISTS (SELECT 1 FROM enumeration e WHERE e.molecule_id = m.id)) TO STDOUT'

    _do_filter(sql, output_file, filters, smarts, dry_run=dry_run, prefix='m.', suffix=suffix)


def filter_need_conf(output_file, count, filters, smarts, dry_run=False):

    utils.expand_path(output_file)

    sql = "COPY (SELECT e.smiles, e.id, e.code FROM enumeration e JOIN molecule m ON e.molecule_id = m.id " + \
        "WHERE NOT EXISTS (SELECT 1 FROM conformer c WHERE c.enumeration_id = e.id) AND"
    if count:
        suffix = ' LIMIT ' + str(int(count)) + ') TO STDOUT'
    else:
        suffix = ') TO STDOUT'

    _do_filter(sql, output_file, filters, smarts, dry_run=dry_run, prefix='m.', suffix=suffix)


def _gen_enumerated_query(filters, smarts, codes=None, count=None):

    filters = _gen_filters(filters, smarts, prefix='m.')

    sql = "SELECT e.id, e.molecule_id, e.smiles e_smiles, e.coords, e.code, m.smiles m_smiles, " + \
          "concat_ws(', ', VARIADIC array_agg(s.code)) codes, " + \
          "concat_ws(', ', VARIADIC array_agg(l.name)) libs " + \
          "FROM enumeration e " + \
          "JOIN molecule m ON m.id = e.molecule_id " + \
          "JOIN supply s on s.molecule_id = m.id " + \
          "JOIN file f ON f.id = s.file_id " + \
          "JOIN library l ON l.id = f.library_id " + \
          "GROUP BY e.id, m.id HAVING" + filters
    if codes:
        sql = sql + " AND e.code IN ('" + "','".join(codes) + "')"
    if count:
        sql = sql + " LIMIT " + str(int(count))

    return sql


def _gen_conformers_query(filters, smarts, codes=None, count=None):

    filters = _gen_filters(filters, smarts, prefix='m.')

    sql = 'SELECT c.id AS c_id, e.id AS e_id, e.molecule_id AS m_id, e.smiles, c.coords, e.code, m.smiles ' +\
          'FROM conformer c ' +\
          'JOIN enumeration e ON c.enumeration_id = e.id JOIN molecule m ON m.id = e.molecule_id WHERE' + filters
    if codes:
        sql = sql + " AND e.code IN ('" + "','".join(codes) + "')"
    if count:
        sql = sql + " LIMIT " + str(int(count))

    return sql


def filter_enumerated_cxsmi(output_file, filters, smarts, count, codes=None, dry_run=False):

    utils.expand_path(output_file)

    sql = _gen_enumerated_query(filters, smarts, codes=codes, count=count)
    utils.log('SQL:\n', sql)

    if not dry_run:
        count = 0
        t0 = time.time()
        with open(output_file, 'wt') as writer:
            with engine.connect() as conn:
                DmLog.emit_event('Executing query')
                results = conn.execute(text(sql))
                for result in results:
                    count += 1
                    # e.id, e.molecule_id, e_smiles, e.coords, e.code, m_smiles, codes, libs
                    #    0              1          2        3       4         5      6     7
                    line = "{} |{}|\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                        .format(result[2], result[3], str(result[0]), str(result[1]), result[4], result[5], result[6], result[7])
                    writer.write(line)
        t1 = time.time()
        DmLog.emit_event('Generated {} records in file {} in {}s'.format(count, output_file, round(t1 - t0)))


def filter_enumerated_sdf(output_file, filters, smarts, count, codes=None, dry_run=False):

    utils.expand_path(output_file)

    sql = _gen_enumerated_query(filters, smarts, codes=codes, count=count)
    utils.log('SQL:', sql)

    if not dry_run:
        count = 0
        t0 = time.time()
        with open(output_file, 'wt') as writer:
            with engine.connect() as conn:
                DmLog.emit_event('Executing query')
                results = conn.execute(text(sql))
                for result in results:
                    count += 1
                    # e.id, e.molecule_id, e_smiles, e.coords, e.code, m_smiles, codes, libs
                    #    0              1          2        3       4         5      6     7
                    mol = Chem.MolFromSmiles(result[2] + ' |' + result[3] + '|')
                    if mol.HasProp('_CXSMILES_Data'):
                        mol.ClearProp('_CXSMILES_Data')
                    mol.SetProp('_Name', str(result[0]))
                    mol.SetProp('enum_smi', result[2])
                    mol.SetProp('std_smi', result[5])
                    mol.SetProp('molecule_id', str(result[1]))
                    mol.SetProp('enum_code', result[4])
                    mol.SetProp('codes', result[6])
                    mol.SetProp('libraries', result[7])

                    sdf_block = Chem.SDWriter.GetText(mol)
                    chg_block = rdkit_utils.updateChargeFlagInAtomBlock(sdf_block)
                    writer.write(chg_block)
        t1 = time.time()
        DmLog.emit_event('Generated {} records in file {} in {}s'.format(count, output_file, round(t1 - t0)))


def filter_conformers_cxsmi(output_file, filters, smarts, count, codes=None, dry_run=False):

    utils.expand_path(output_file)

    sql = _gen_conformers_query(filters, smarts, codes=codes, count=count)
    utils.log('SQL:', sql)

    if not dry_run:
        count = 0
        t0 = time.time()
        with open(output_file, 'wt') as writer:
            with engine.connect() as conn:
                DmLog.emit_event('Executing query')
                results = conn.execute(text(sql))
                for result in results:
                    count += 1
                    # SELECT c.id AS c_id, e.id AS e_id, e.molecule_id AS m_id, e.smiles, c.coords, e.code, m.smiles
                    #                0             1                      2     3         4         5       6
                    line = "{} {}\t{}-{}-{}\t{}\t{}\n".format(
                        result[3], result[4], str(result[2]), str(result[1]), str(result[0]), result[6], result[5])
                    writer.write(line)
        t1 = time.time()
        DmLog.emit_event('Generated {} records in file {} in {}s'.format(count, output_file, round(t1 - t0)))


def filter_conformers_sdf(output_file, filters, smarts, count, codes=None, dry_run=False):

    utils.expand_path(output_file)

    sql = _gen_conformers_query(filters, smarts, codes=codes, count=count)
    utils.log('SQL:', sql)

    if not dry_run:
        count = 0
        t0 = time.time()
        with open(output_file, 'wt') as writer:
            with engine.connect() as conn:
                DmLog.emit_event('Executing query')
                results = conn.execute(text(sql))
                # SELECT c.id AS c_id, e.id AS e_id, e.molecule_id AS m_id, e.smiles, c.coords, e.code, m.smiles
                #                0             1                      2     3         4         5       6

                for result in results:
                    count += 1
                    mol = Chem.MolFromSmiles(result[3] + ' ' + result[4])
                    if mol.HasProp('_CXSMILES_Data'):
                        mol.ClearProp('_CXSMILES_Data')
                    id = '{}-{}-{}'.format(result[2], result[1], result[0])
                    mol.SetProp('_Name', id)
                    mol.SetProp('enum_smi', result[3])
                    mol.SetProp('std_smi', result[6])
                    mol.SetProp('enum_code', result[5])
                    mol.SetIntProp('m_id', result[2])
                    mol.SetIntProp('e_id', result[1])
                    mol.SetIntProp('c_id', result[0])

                    sdf_block = Chem.SDWriter.GetText(mol)
                    chg_block = rdkit_utils.updateChargeFlagInAtomBlock(sdf_block)
                    writer.write(chg_block)
        t1 = time.time()
        DmLog.emit_event('Generated {} records in file {} in {}s'.format(count, output_file, round(t1 - t0)))


def main():

    # Example:
    #   python -m moldb.filter --output-molecules filtered.smi --output-need-enum need-enum.smi --specification specification.txt

    parser = argparse.ArgumentParser(description='Filter')
    parser.add_argument('-s', '--specification', help="Filter specification file")
    parser.add_argument('--output-molecules', help="Output file for SMILES")
    parser.add_argument('--output-need-enum', help="Output file for SMILES needing enumeration")
    parser.add_argument('--output-need-conf', help="Output file for SMILES needing conformer generation")
    parser.add_argument('--output-enumerated', help="Output SDF or CXSMILES file for enumerated molecules")
    parser.add_argument('--output-conformer', help="Output SDF file for 3D conformers")
    parser.add_argument('-m', '--enum-codes', nargs='+', help="Enumerated types (B, C, T, M)")
    parser.add_argument('-c', '--count', type=int, help='Max number of molecules to extract')
    parser.add_argument('--dry-run', action='store_true', help="Generate SQL but don't execute")

    args = parser.parse_args()
    DmLog.emit_event("filter: ", args)

    if args.specification:
        filters, smarts = moldb_utils.read_specification(args.specification)
    else:
        filters = {}
        smarts = []

    if args.output_molecules:
        filter_molecules(args.output_molecules, args.count, filters, smarts, dry_run=args.dry_run)

    if args.output_need_enum:
        filter_need_enum(args.output_need_enum, args.count, filters, smarts, dry_run=args.dry_run)

    if args.output_need_conf:
        filter_need_conf(args.output_need_conf, args.count, filters, smarts, dry_run=args.dry_run)

    if args.output_enumerated:
        if args.output_enumerated.endswith('.cxsmi'):
            filter_enumerated_cxsmi(args.output_enumerated, filters, smarts, args.count, codes=args.enum_codes, dry_run=args.dry_run)
        elif args.output_enumerated.endswith('.sdf'):
            filter_enumerated_sdf(args.output_enumerated, filters, smarts, args.count, codes=args.enum_codes, dry_run=args.dry_run)
        else:
            raise ValueError('Output file must be .sdf or .cxsmi')

    if args.output_conformer:
        if args.output_conformer.endswith('.cxsmi'):
            filter_conformers_cxsmi(args.output_conformer, filters, smarts, args.count, codes=args.enum_codes, dry_run=args.dry_run)
        elif args.output_conformer.endswith('.sdf'):
            filter_conformers_sdf(args.output_conformer, filters, smarts, args.count, codes=args.enum_codes, dry_run=args.dry_run)
        else:
            raise ValueError('Output file must be .sdf or .cxsmi')


if __name__ == "__main__":
    main()
