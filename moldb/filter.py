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


def _gen_filters(
            min_hac=None, max_hac=None,
            min_rotb=None, max_rotb=None,
            min_rings=None, max_rings=None,
            min_aro_rings=None, max_aro_rings=None,
            min_chiral_centres=None, max_chiral_centres=None,
            min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
            min_sp3=None, max_sp3=None,
            min_logp=None, max_logp=None,
            min_tpsa=None, max_tpsa=None,
            prefix=''):

    terms = (
        (min_hac, max_hac, 'hac'),
        (min_rotb, max_rotb, 'rot_bonds'),
        (min_rings, max_rings, 'ring_count'),
        (min_aro_rings, max_aro_rings, 'aromatic_ring_count'),
        (min_chiral_centres, max_chiral_centres, 'chiral_centres'),
        (min_undefined_chiral_centres, max_undefined_chiral_centres, 'undef_chiral_centres'),
        (min_sp3, max_sp3, 'num_sp3'),
        (min_logp, max_logp, 'logp'),
        (min_tpsa, max_tpsa, 'tpsa')
    )
    sql = " "
    count = 0
    for term in terms:
        sql, count = _add_term(sql, term[0], term[1], term[2], count, prefix=prefix)

    return sql


def _do_filter(sql, output_file, dry_run=False, prefix='', suffix='',
               min_hac=None, max_hac=None,
               min_rotb=None, max_rotb=None,
               min_rings=None, max_rings=None,
               min_aro_rings=None, max_aro_rings=None,
               min_chiral_centres=None, max_chiral_centres=None,
               min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
               min_sp3=None, max_sp3=None,
               min_logp=None, max_logp=None,
               min_tpsa=None, max_tpsa=None):

    utils.expand_path(output_file)

    filters = _gen_filters(prefix=prefix,
                           min_hac=min_hac, max_hac=max_hac,
                           min_rotb=min_rotb, max_rotb=max_rotb,
                           min_rings=min_rings, max_rings=max_rings,
                           min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
                           min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
                           min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
                           min_sp3=min_sp3, max_sp3=max_sp3,
                           min_logp=min_logp, max_logp=max_logp,
                           min_tpsa=min_tpsa, max_tpsa=max_tpsa)

    sql = sql + filters + suffix

    DmLog.emit_event('Query SQL:\n', sql)

    if not dry_run:
        with Session(engine) as session:

            output = open(output_file, 'w')
            cursor = session.connection().connection.cursor()
            cursor.copy_expert(sql, output)

            session.commit()


def filter_smiles(output_file, count, dry_run=False,
           min_hac=None, max_hac=None,
           min_rotb=None, max_rotb=None,
           min_rings=None, max_rings=None,
           min_aro_rings=None, max_aro_rings=None,
           min_chiral_centres=None, max_chiral_centres=None,
           min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
           min_sp3=None, max_sp3=None,
           min_logp=None, max_logp=None,
           min_tpsa=None, max_tpsa=None):


    sql = "COPY (SELECT smiles, id FROM molecule WHERE "
    if count:
        suffix = ' LIMIT ' + str(int(count)) + ') TO STDOUT'
    else:
        suffix = ') TO STDOUT'

    _do_filter(sql, output_file, dry_run=dry_run, suffix=suffix,
               min_hac=min_hac, max_hac=max_hac,
               min_rotb=min_rotb, max_rotb=max_rotb,
               min_rings=min_rings, max_rings=max_rings,
               min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
               min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
               min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
               min_sp3=min_sp3, max_sp3=max_sp3,
               min_logp=min_logp, max_logp=max_logp,
               min_tpsa=min_tpsa, max_tpsa=max_tpsa)


def filter_need_enum(output_file, count, dry_run=False,
                     min_hac=None, max_hac=None,
                     min_rotb=None, max_rotb=None,
                     min_rings=None, max_rings=None,
                     min_aro_rings=None, max_aro_rings=None,
                     min_chiral_centres=None, max_chiral_centres=None,
                     min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
                     min_sp3=None, max_sp3=None,
                     min_logp=None, max_logp=None,
                     min_tpsa=None, max_tpsa=None):

    utils.expand_path(output_file)

    sql = "COPY (SELECT smiles, id FROM molecule m WHERE"
    if count:
        suffix = ' AND id NOT IN (SELECT DISTINCT molecule_id FROM enumeration) LIMIT ' + str(int(count)) + ') TO STDOUT'
    else:
        suffix = ' AND id NOT IN (SELECT DISTINCT molecule_id FROM enumeration)) TO STDOUT'

    _do_filter(sql, output_file, dry_run=dry_run, suffix=suffix,
               min_hac=min_hac, max_hac=max_hac,
               min_rotb=min_rotb, max_rotb=max_rotb,
               min_rings=min_rings, max_rings=max_rings,
               min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
               min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
               min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
               min_sp3=min_sp3, max_sp3=max_sp3,
               min_logp=min_logp, max_logp=max_logp,
               min_tpsa=min_tpsa, max_tpsa=max_tpsa)


def filter_need_conf(output_file, count, dry_run=False,
                     min_hac=None, max_hac=None,
                     min_rotb=None, max_rotb=None,
                     min_rings=None, max_rings=None,
                     min_aro_rings=None, max_aro_rings=None,
                     min_chiral_centres=None, max_chiral_centres=None,
                     min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
                     min_sp3=None, max_sp3=None,
                     min_logp=None, max_logp=None,
                     min_tpsa=None, max_tpsa=None):

    utils.expand_path(output_file)

    sql = "COPY (SELECT e.smiles, e.id, e.code FROM enumeration e " + \
        "WHERE e.molecule_id NOT IN " + \
          "(SELECT DISTINCT e.id FROM enumeration e JOIN molecule m ON e.molecule_id = m.id WHERE"
    if count:
        suffix = ') LIMIT ' + str(int(count)) + ') TO STDOUT'
    else:
        suffix = ')) TO STDOUT'

    _do_filter(sql, output_file, dry_run=dry_run, prefix='m.', suffix=suffix,
               min_hac=min_hac, max_hac=max_hac,
               min_rotb=min_rotb, max_rotb=max_rotb,
               min_rings=min_rings, max_rings=max_rings,
               min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
               min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
               min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
               min_sp3=min_sp3, max_sp3=max_sp3,
               min_logp=min_logp, max_logp=max_logp,
               min_tpsa=min_tpsa, max_tpsa=max_tpsa
               )


def _gen_enumerated_query(codes=None, count=None,
                          min_hac=None, max_hac=None,
                          min_rotb=None, max_rotb=None,
                          min_rings=None, max_rings=None,
                          min_aro_rings=None, max_aro_rings=None,
                          min_chiral_centres=None, max_chiral_centres=None,
                          min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
                          min_sp3=None, max_sp3=None,
                          min_logp=None, max_logp=None,
                          min_tpsa=None, max_tpsa=None):

    filters = _gen_filters(min_hac=min_hac, max_hac=max_hac,
                           min_rotb=min_rotb, max_rotb=max_rotb,
                           min_rings=min_rings, max_rings=max_rings,
                           min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
                           min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
                           min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
                           min_sp3=min_sp3, max_sp3=max_sp3,
                           min_logp=min_logp, max_logp=max_logp,
                           min_tpsa=min_tpsa, max_tpsa=max_tpsa,
                           prefix='m.')

    sql = 'SELECT e.molecule_id, e.smiles, e.code, m.smiles FROM enumeration e JOIN molecule m ON m.id = e.molecule_id WHERE' + filters
    if codes:
        sql = sql + " AND e.code IN ('" + "','".join(codes) + "')"
    if count:
        sql = sql + " LIMIT " + str(int(count))

    DmLog.emit_event('Query SQL:\n', sql)
    return sql


def _gen_conformers_query(codes=None, count=None,
                          min_hac=None, max_hac=None,
                          min_rotb=None, max_rotb=None,
                          min_rings=None, max_rings=None,
                          min_aro_rings=None, max_aro_rings=None,
                          min_chiral_centres=None, max_chiral_centres=None,
                          min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
                          min_sp3=None, max_sp3=None,
                          min_logp=None, max_logp=None,
                          min_tpsa=None, max_tpsa=None):

    filters = _gen_filters(min_hac=min_hac, max_hac=max_hac,
                           min_rotb=min_rotb, max_rotb=max_rotb,
                           min_rings=min_rings, max_rings=max_rings,
                           min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
                           min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
                           min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
                           min_sp3=min_sp3, max_sp3=max_sp3,
                           min_logp=min_logp, max_logp=max_logp,
                           min_tpsa=min_tpsa, max_tpsa=max_tpsa,
                           prefix='m.')

    sql = 'SELECT c.id AS c_id, e.id AS e_id, e.molecule_id AS m_id, e.smiles, c.coords, e.code, m.smiles FROM conformer c ' \
          'JOIN enumeration e ON c.enumeration_id = e.id JOIN molecule m ON m.id = e.molecule_id WHERE' + filters
    if codes:
        sql = sql + " AND e.code IN ('" + "','".join(codes) + "')"
    if count:
        sql = sql + " LIMIT " + str(int(count))

    DmLog.emit_event('Query SQL:\n', sql)
    return sql


def filter_enumerated_cxsmi(output_file, count, codes=None,
                     min_hac=None, max_hac=None,
                     min_rotb=None, max_rotb=None,
                     min_rings=None, max_rings=None,
                     min_aro_rings=None, max_aro_rings=None,
                     min_chiral_centres=None, max_chiral_centres=None,
                     min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
                     min_sp3=None, max_sp3=None,
                     min_logp=None, max_logp=None,
                     min_tpsa=None, max_tpsa=None):

    utils.expand_path(output_file)

    sql = _gen_enumerated_query(codes=codes, count=count,
                                min_hac=min_hac, max_hac=max_hac,
                                min_rotb=min_rotb, max_rotb=max_rotb,
                                min_rings=min_rings, max_rings=max_rings,
                                min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
                                min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
                                min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
                                min_sp3=min_sp3, max_sp3=max_sp3,
                                min_logp=min_logp, max_logp=max_logp,
                                min_tpsa=min_tpsa, max_tpsa=max_tpsa)

    count = 0
    t0 = time.time()
    with open(output_file, 'wt') as writer:
        with engine.connect() as conn:
            results = conn.execute(text(sql))
            for result in results:
                count += 1
                line = "{}\t{}\t{}\t{}\n".format(result[1], str(result[0]), result[2], result[3])
                writer.write(line)
    t1 = time.time()
    DmLog.emit_event('Generated {} records in file {} in {}s'.format(count, output_file, round(t1 - t0)))


def filter_enumerated_sdf(output_file, count, codes=None,
                            min_hac=None, max_hac=None,
                            min_rotb=None, max_rotb=None,
                            min_rings=None, max_rings=None,
                            min_aro_rings=None, max_aro_rings=None,
                            min_chiral_centres=None, max_chiral_centres=None,
                            min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
                            min_sp3=None, max_sp3=None,
                            min_logp=None, max_logp=None,
                            min_tpsa=None, max_tpsa=None):

    utils.expand_path(output_file)

    sql = _gen_enumerated_query(codes=codes, count=count,
                                min_hac=min_hac, max_hac=max_hac,
                                min_rotb=min_rotb, max_rotb=max_rotb,
                                min_rings=min_rings, max_rings=max_rings,
                                min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
                                min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
                                min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
                                min_sp3=min_sp3, max_sp3=max_sp3,
                                min_logp=min_logp, max_logp=max_logp,
                                min_tpsa=min_tpsa, max_tpsa=max_tpsa)

    count = 0
    t0 = time.time()
    with open(output_file, 'wt') as writer:
        with engine.connect() as conn:
            results = conn.execute(text(sql))
            for result in results:
                count += 1
                mol = Chem.MolFromSmiles(result[1])
                enum_smi = result[1].split(' ')[0]
                if mol.HasProp('_CXSMILES_Data'):
                    mol.ClearProp('_CXSMILES_Data')
                mol.SetProp('_Name', str(result[0]))
                mol.SetProp('enum_smi', enum_smi)
                mol.SetProp('std_smi', result[3])
                mol.SetProp('enum_code', result[2])

                sdf_block = Chem.SDWriter.GetText(mol)
                chg_block = rdkit_utils.updateChargeFlagInAtomBlock(sdf_block)
                writer.write(chg_block)
    t1 = time.time()
    DmLog.emit_event('Generated {} records in file {} in {}s'.format(count, output_file, round(t1 - t0)))


def filter_conformers_cxsmi(output_file, count, codes=None,
                            min_hac=None, max_hac=None,
                            min_rotb=None, max_rotb=None,
                            min_rings=None, max_rings=None,
                            min_aro_rings=None, max_aro_rings=None,
                            min_chiral_centres=None, max_chiral_centres=None,
                            min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
                            min_sp3=None, max_sp3=None,
                            min_logp=None, max_logp=None,
                            min_tpsa=None, max_tpsa=None):

    utils.expand_path(output_file)

    sql = _gen_conformers_query(codes=codes, count=count,
                                min_hac=min_hac, max_hac=max_hac,
                                min_rotb=min_rotb, max_rotb=max_rotb,
                                min_rings=min_rings, max_rings=max_rings,
                                min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
                                min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
                                min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
                                min_sp3=min_sp3, max_sp3=max_sp3,
                                min_logp=min_logp, max_logp=max_logp,
                                min_tpsa=min_tpsa, max_tpsa=max_tpsa)

    count = 0
    t0 = time.time()
    with open(output_file, 'wt') as writer:
        with engine.connect() as conn:
            results = conn.execute(text(sql))
            for result in results:
                count += 1
                # SELECT c.id AS c_id, e.id AS e_id, e.molecule_id AS m_id, e.smiles, c.coords, e.code, m.smiles
                #                0             1                      2     3         4         5       6
                line = "{} {}\t{}-{}-{}\t{}\t{}\n".format(result[3], result[4], str(result[2]), str(result[1]), str(result[0]), result[6], result[5])
                writer.write(line)
    t1 = time.time()
    DmLog.emit_event('Generated {} records in file {} in {}s'.format(count, output_file, round(t1 - t0)))


def filter_conformers_sdf(output_file, count, codes=None,
                          min_hac=None, max_hac=None,
                          min_rotb=None, max_rotb=None,
                          min_rings=None, max_rings=None,
                          min_aro_rings=None, max_aro_rings=None,
                          min_chiral_centres=None, max_chiral_centres=None,
                          min_undefined_chiral_centres=None, max_undefined_chiral_centres=None,
                          min_sp3=None, max_sp3=None,
                          min_logp=None, max_logp=None,
                          min_tpsa=None, max_tpsa=None):

    utils.expand_path(output_file)

    sql = _gen_conformers_query(codes=codes, count=count,
                                min_hac=min_hac, max_hac=max_hac,
                                min_rotb=min_rotb, max_rotb=max_rotb,
                                min_rings=min_rings, max_rings=max_rings,
                                min_aro_rings=min_aro_rings, max_aro_rings=max_aro_rings,
                                min_chiral_centres=min_chiral_centres, max_chiral_centres=max_chiral_centres,
                                min_undefined_chiral_centres=min_undefined_chiral_centres, max_undefined_chiral_centres=max_undefined_chiral_centres,
                                min_sp3=min_sp3, max_sp3=max_sp3,
                                min_logp=min_logp, max_logp=max_logp,
                                min_tpsa=min_tpsa, max_tpsa=max_tpsa)

    count = 0
    t0 = time.time()
    with open(output_file, 'wt') as writer:
        with engine.connect() as conn:
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
    #   python -m moldb.filter --output-smiles filtered.smi --output-need-enum need-enum.smi --min-hac 16 --max-hac 24 \
    #     --min-rings 2 --min-aro-rings 1 --max-chiral-centres 2 --max-undefined-chiral-centres 0

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Filter')
    parser.add_argument('--output-smiles', help="Output file for SMILES")
    parser.add_argument('--output-need-enum', help="Output file for SMILES needing enumeration")
    parser.add_argument('--output-need-conf', help="Output file for SMILES needing conformer generation")
    parser.add_argument('--output-enumerated', help="Output SDF or CXSMILES file for enumerated molecules")
    parser.add_argument('--output-conformer', help="Output SDF file for 3D conformers")

    moldb_utils.add_filter_args(parser)

    parser.add_argument('-m', '--enum-codes', nargs='+', help="Enumerated types (B, C, T, M)")
    parser.add_argument('-c', '--count', type=int, help='Max number of molecules to extract')
    parser.add_argument('--dry-run', action='store_true', help="Generate SQL but don't execute")

    args = parser.parse_args()
    DmLog.emit_event("filter: ", args)

    if args.output_smiles:
        filter_smiles(args.output_smiles, args.count, dry_run=args.dry_run,
            min_hac=args.min_hac, max_hac=args.max_hac,
            min_rotb=args.min_rotb, max_rotb=args.max_rotb,
            min_rings=args.min_rings, max_rings=args.max_rings,
            min_aro_rings=args.min_aro_rings, max_aro_rings=args.max_aro_rings,
            min_chiral_centres=args.min_chiral_centres, max_chiral_centres=args.max_chiral_centres,
            min_undefined_chiral_centres=args.min_undefined_chiral_centres,
                      max_undefined_chiral_centres=args.max_undefined_chiral_centres,
            min_sp3=args.min_sp3, max_sp3=args.max_sp3,
            min_logp=args.min_logp, max_logp=args.max_logp,
            min_tpsa=args.min_tpsa, max_tpsa=args.max_tpsa)

    if args.output_need_enum:
        filter_need_enum(args.output_need_enum, args.count, dry_run=args.dry_run,
                          min_hac=args.min_hac, max_hac=args.max_hac,
                          min_rotb=args.min_rotb, max_rotb=args.max_rotb,
                          min_rings=args.min_rings, max_rings=args.max_rings,
                          min_aro_rings=args.min_aro_rings, max_aro_rings=args.max_aro_rings,
                          min_chiral_centres=args.min_chiral_centres, max_chiral_centres=args.max_chiral_centres,
                          min_undefined_chiral_centres=args.min_undefined_chiral_centres,
                         max_undefined_chiral_centres=args.max_undefined_chiral_centres,
                          min_sp3=args.min_sp3, max_sp3=args.max_sp3,
                          min_logp=args.min_logp, max_logp=args.max_logp,
                          min_tpsa=args.min_tpsa, max_tpsa=args.max_tpsa)

    if args.output_need_conf:
        filter_need_conf(args.output_need_conf, args.count, dry_run=args.dry_run,
                         min_hac=args.min_hac, max_hac=args.max_hac,
                         min_rotb=args.min_rotb, max_rotb=args.max_rotb,
                         min_rings=args.min_rings, max_rings=args.max_rings,
                         min_aro_rings=args.min_aro_rings, max_aro_rings=args.max_aro_rings,
                         min_chiral_centres=args.min_chiral_centres, max_chiral_centres=args.max_chiral_centres,
                         min_undefined_chiral_centres=args.min_undefined_chiral_centres,
                         max_undefined_chiral_centres=args.max_undefined_chiral_centres,
                         min_sp3=args.min_sp3, max_sp3=args.max_sp3,
                         min_logp=args.min_logp, max_logp=args.max_logp,
                         min_tpsa=args.min_tpsa, max_tpsa=args.max_tpsa)

    if args.output_enumerated:
        if args.output_enumerated.endswith('.cxsmi'):
            filter_enumerated_cxsmi(args.output_enumerated, args.count, codes=args.enum_codes,
                                min_hac=args.min_hac, max_hac=args.max_hac,
                                min_rotb=args.min_rotb, max_rotb=args.max_rotb,
                                min_rings=args.min_rings, max_rings=args.max_rings,
                                min_aro_rings=args.min_aro_rings, max_aro_rings=args.max_aro_rings,
                                min_chiral_centres=args.min_chiral_centres, max_chiral_centres=args.max_chiral_centres,
                                min_undefined_chiral_centres=args.min_undefined_chiral_centres,
                                max_undefined_chiral_centres=args.max_undefined_chiral_centres,
                                min_sp3=args.min_sp3, max_sp3=args.max_sp3,
                                min_logp=args.min_logp, max_logp=args.max_logp,
                                min_tpsa=args.min_tpsa, max_tpsa=args.max_tpsa)
        elif args.output_enumerated.endswith('.sdf'):
            filter_enumerated_sdf(args.output_enumerated, args.count, codes=args.enum_codes,
                              min_hac=args.min_hac, max_hac=args.max_hac,
                              min_rotb=args.min_rotb, max_rotb=args.max_rotb,
                              min_rings=args.min_rings, max_rings=args.max_rings,
                              min_aro_rings=args.min_aro_rings, max_aro_rings=args.max_aro_rings,
                              min_chiral_centres=args.min_chiral_centres, max_chiral_centres=args.max_chiral_centres,
                              min_undefined_chiral_centres=args.min_undefined_chiral_centres,
                              max_undefined_chiral_centres=args.max_undefined_chiral_centres,
                              min_sp3=args.min_sp3, max_sp3=args.max_sp3,
                              min_logp=args.min_logp, max_logp=args.max_logp,
                              min_tpsa=args.min_tpsa, max_tpsa=args.max_tpsa)
        else:
            raise ValueError('Output file must be .sdf or .cxsmi')

    if args.output_conformer:
        if args.output_conformer.endswith('.cxsmi'):
            filter_conformers_cxsmi(args.output_conformer, args.count, codes=args.enum_codes,
                                    min_hac=args.min_hac, max_hac=args.max_hac,
                                    min_rotb=args.min_rotb, max_rotb=args.max_rotb,
                                    min_rings=args.min_rings, max_rings=args.max_rings,
                                    min_aro_rings=args.min_aro_rings, max_aro_rings=args.max_aro_rings,
                                    min_chiral_centres=args.min_chiral_centres, max_chiral_centres=args.max_chiral_centres,
                                    min_undefined_chiral_centres=args.min_undefined_chiral_centres,
                                    max_undefined_chiral_centres=args.max_undefined_chiral_centres,
                                    min_sp3=args.min_sp3, max_sp3=args.max_sp3,
                                    min_logp=args.min_logp, max_logp=args.max_logp,
                                    min_tpsa=args.min_tpsa, max_tpsa=args.max_tpsa)
        elif args.output_conformer.endswith('.sdf'):
            filter_conformers_sdf(args.output_conformer, args.count, codes=args.enum_codes,
                                  min_hac=args.min_hac, max_hac=args.max_hac,
                                  min_rotb=args.min_rotb, max_rotb=args.max_rotb,
                                  min_rings=args.min_rings, max_rings=args.max_rings,
                                  min_aro_rings=args.min_aro_rings, max_aro_rings=args.max_aro_rings,
                                  min_chiral_centres=args.min_chiral_centres, max_chiral_centres=args.max_chiral_centres,
                                  min_undefined_chiral_centres=args.min_undefined_chiral_centres,
                                  max_undefined_chiral_centres=args.max_undefined_chiral_centres,
                                  min_sp3=args.min_sp3, max_sp3=args.max_sp3,
                                  min_logp=args.min_logp, max_logp=args.max_logp,
                                  min_tpsa=args.min_tpsa, max_tpsa=args.max_tpsa)
        else:
            raise ValueError('Output file must be .sdf or .cxsmi')

if __name__ == "__main__":
    main()
