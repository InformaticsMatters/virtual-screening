#!/usr/bin/env python

# Copyright 2024 Informatics Matters Ltd.
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

import argparse, datetime, time, sys, os

from mordred import Calculator, descriptors

from rdkit import Chem

import utils, rdkit_utils
from dm_job_utilities.dm_log import DmLog


def run(input, outfile, mode='hac', include_3d=False, delimiter=None, id_column=None, omit_fields=False,
        read_header=False, write_header=False,
        read_records=50, interval=1000):

    if include_3d:
        if not (input.endswith('.sdf') or input.endswith('.sdf.gz')):
            DmLog.emit_event('If calculating 3D descriptors input must be a SDF with 3D molecules')
            exit(1)
        if not rdkit_utils.check_molecules_are_3d(input):
            DmLog.emit_event('Molecules do not seem to have 3D coordinates')
            exit(1)

    utils.expand_path(outfile)

    calc = Calculator(descriptors, ignore_3D=not include_3d)

    count = 0
    errors = 0

    t0 = time.time()

    # setup the reader
    reader = rdkit_utils.create_reader(input, id_column=id_column, read_records=read_records,
                                       read_header=read_header, delimiter=delimiter)
    extra_field_names = reader.get_extra_field_names()

    calc_field_names = [str(name) for name in calc.descriptors]
    DmLog.emit_event("Calculating {} descriptors".format(len(calc_field_names)))

    # setup the writer
    writer = rdkit_utils.create_writer(outfile,
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
            DmLog.emit_event("Failed to read molecule for record", count)
            continue

        if omit_fields:
            for name in mol.GetPropNames():
                mol.ClearProp(name)
            props = []

        try:
            if mode == 'none':
                biggest = mol
                cann_smi = smi
            else:
                biggest = rdkit_utils.fragment(mol, mode)
                cann_smi = Chem.MolToSmiles(biggest)

            values = calc(biggest)

        except KeyboardInterrupt:
            utils.log('Interrupted')
            sys.exit(0)

        except:
            errors += 1
            DmLog.emit_event('Failed to process record', count)
            continue

        writer.write(cann_smi, biggest, id, props, values)

    t1 = time.time()
    utils.log('Processing took {} secs'.format(round(t1 - t0)))


def main():

    # Example:
    #   python -m im_mordred.descriptor_generator -i data/10.smi -o descriptors.smi -d tab

    # ----- command line args definitions ---------------------------------------------

    parser = argparse.ArgumentParser(description='Mordred 2D descriptors')
    parser.add_argument('-i', '--input', required=True, help="Input file (.smi or .sdf)")
    parser.add_argument('-o', '--output', default='descriptors2d.smi', help="Output file (.smi or .sdf")

    parser.add_argument('-f', '--fragment-method', choices=['hac', 'mw', 'none'], default='hac',
                        help='Strategy for picking largest fragment (mw or hac or none')

    parser.add_argument('-3', '--include-3d', action='store_true',
                        help="Include 3D descriptors (requires 3D molecules in SDF file)")

    parser.add_argument('-k', '--omit-fields', action='store_true',
                        help="Don't include fields from the input in the output")

    # to pass tab as the delimiter specify it as $'\t' or use one of the symbolic names 'comma', 'tab', 'space' or 'pipe'
    parser.add_argument('-d', '--delimiter', help="Delimiter when using SMILES")
    parser.add_argument('--id-column', help="Column for name field (zero based integer for .smi, text for SDF)")
    parser.add_argument('--read-header', action='store_true',
                        help="Read a header line with the field names when reading .smi or .txt")
    parser.add_argument('--write-header', action='store_true', help='Write a header line when writing .smi or .txt')
    parser.add_argument('--read-records', default=100, type=int,
                        help="Read this many records to determine the fields that are present")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("descriptor_calc: ", args)

    delimiter = utils.read_delimiter(args.delimiter)

    run(args.input, args.output, mode=args.fragment_method, include_3d=args.include_3d,
        omit_fields=args.omit_fields ,delimiter=delimiter, id_column=args.id_column,
        read_header=args.read_header, write_header=args.write_header,
        read_records=args.read_records, interval=args.interval)


if __name__ == "__main__":
    main()
