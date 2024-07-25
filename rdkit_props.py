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


import argparse, os, gzip, time

import utils, rdkit_utils, rdkit_calcs

from dm_job_utilities.dm_log import DmLog


calc_props = {
    # calc key: (result prop name, prop description, function to call, index in result if multiple results)
    'hac': ('RDK_hac', 'Calculate heavy atom count', rdkit_calcs.calc_hac, None),
    'num_rot_bonds': ('RDK_numRotBonds', 'Calculate number of rotatable bonds', rdkit_calcs.calc_num_rot_bonds, None),
    'num_rings': ('RDK_numRings', 'Calculate number of rings', rdkit_calcs.calc_num_rings, None),
    'num_aro_rings': ('RDK_numAroRings', 'Calculate number of aromatic rings', rdkit_calcs.calc_num_aro_rings, None),
    'num_cc': ('RDK_numChiralCentres', 'Calculate number of chiral centres', rdkit_calcs.calc_num_cc, 0),
    'num_undef_cc': ('RDK_numUndefChiralCentres', 'Calculate number of undefined chiral centres', rdkit_calcs.calc_num_cc, 1),
    'num_sp3': ('RDK_numSP3', 'Calculate number of sp3 hybridised carbon atoms', rdkit_calcs.calc_num_sp3, None),
    'logp': ('RDK_logp', 'Calculate logP', rdkit_calcs.calc_logp_crippen, None),
    'tpsa': ('RDK_tpsa', 'Calculate topological polar surface area', rdkit_calcs.calc_tpsa, None),
    'sa_score': ('SAScore', 'Synthetic accessibility score', rdkit_calcs.calc_sa_score, None)
}


def process(input,
            outfile,
            calcs,
            delimiter,
            id_column=None,
            mol_column=0,
            omit_fields=False,
            read_header=False,
            write_header=False,
            read_records=100,
            interval=0):

    utils.log('Using calculations:', calcs)
    for calc in calcs:
        if calc not in calc_props:
            utils.log('calculation', calc, 'is not supported')

    utils.expand_path(outfile)

    count = 0
    errors = 0

    # setup the reader
    reader = rdkit_utils.create_reader(input,
                                       id_column=id_column,
                                       mol_column=mol_column,
                                       read_records=read_records,
                                       read_header=read_header,
                                       delimiter=delimiter)
    extra_field_names = reader.get_extra_field_names()

    calc_field_names = [calc_props[s][0] for s in calcs]

    # setup the writer
    writer = rdkit_utils.create_writer(outfile,
                                       id_column=id_column,
                                       mol_column=mol_column,
                                       extra_field_names=extra_field_names,
                                       calc_prop_names=calc_field_names,
                                       delimiter=delimiter)

    id_col_type, id_col_value = utils.is_type(id_column, int)

    # read the input records and write the output
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break

        mol, smi, id, props = t
        if count == 0 and write_header:
            headers = rdkit_utils.generate_headers(
                id_col_type,
                id_col_value,
                reader.get_mol_field_name(),
                reader.field_names,
                calc_field_names,
                omit_fields)

            writer.write_header(headers)

        count += 1

        if interval and count % interval == 0:
            DmLog.emit_event("Processed {} records".format(count))
        if count % 50000 == 0:
            # Emit a 'total' cost, replacing all prior costs
            DmLog.emit_cost(count)

        if not mol:
            errors += 1
            DmLog.emit_event("Failed to process record", count)
            continue

        # calculate the molecular props
        try:
            values = []
            results_dict = {}
            for calc in calcs:
                if calc in calc_props:
                    tup = calc_props[calc]
                    func = tup[2]
                    result = results_dict.get(func)
                    if not result:
                        result = func(mol)
                        results_dict[func] = result
                    if tup[3] is None:
                        values.append(result)
                    else:
                        values.append(result[tup[3]])

        except Exception as e:
            errors += 1
            DmLog.emit_event('Failed to process record', count, str(e))
            continue

        if omit_fields:
            for name in mol.GetPropNames():
                mol.ClearProp(name)
            props = []

        # write the output
        writer.write(smi, mol, id, props, values)

    writer.close()
    reader.close()

    return count, errors


def main():
    # Example usage:
    #   ./rdkit_props.py -i data/1000.smi --outfile out.sdf -a --delimiter tab --id-column 1 --interval 100
    #   ./rdkit_props.py -i data/1000.smi --outfile out.smi --write-header -c logp tpsa --delimiter tab --id-column 1 --interval 100

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Calc RDKit props')
    parser.add_argument('-i', '--input', required=True, help="Input file as SMILES or SDF")
    parser.add_argument('-o', '--outfile', required=True, help="Output file as SMILES or SDF")

    parser.add_argument('-k', '--omit-fields', action='store_true',
                        help="Don't include fields from the input in the output")

    # to pass tab as the delimiter specify it as $'\t' or use one of the symbolic names 'comma', 'tab', 'space' or 'pipe'
    parser.add_argument('-d', '--delimiter', help="Delimiter when using SMILES")
    parser.add_argument('--id-column', help="Column for name field (zero based integer for .smi, text for SDF)")
    parser.add_argument('--mol-column', type=int, default=0,
                        help="Column index for molecule when using delineated text formats (zero based integer)")
    parser.add_argument('--read-header', action='store_true',
                        help="Read a header line with the field names when reading .smi or .txt")
    parser.add_argument('--write-header', action='store_true', help='Write a header line when writing .smi or .txt')
    parser.add_argument('--read-records', default=100, type=int,
                        help="Read this many records to determine the fields that are present")

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-a', '--all', action='store_true', help="Calculate all properties")
    calc_choices = [k for k in calc_props]
    group.add_argument('-c', '--calcs', nargs='+', choices=calc_choices, help="Calculators to use")

    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("rdk_props.py: ", args)

    if args.all:
        calcs_to_use = calc_props.keys()
    else:
        calcs_to_use = args.calcs

    # special processing of delimiter to allow it to be set as a name

    delimiter = utils.read_delimiter(args.delimiter)

    t0 = time.time()
    count, errors = process(
        args.input, args.outfile, calcs_to_use, delimiter, id_column=args.id_column, mol_column=args.mol_column,
        omit_fields=args.omit_fields, read_header=args.read_header, write_header=args.write_header,
        read_records=args.read_records, interval=args.interval, )
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processed {} records in {} seconds. {} errors.'.format(count, duration_s, errors))
    # Emit final 'total' cost, replacing all prior costs
    DmLog.emit_cost(count * len(calcs_to_use))


if __name__ == "__main__":
    main()
