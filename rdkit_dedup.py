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


import argparse, os, gzip, time
import traceback

import utils, rdkit_utils
from dm_job_utilities.dm_log import DmLog

from rdkit import Chem


def process(input,
            outfile,
            mode='hac',
            delimiter=None,
            id_column=None,
            mol_column=0,
            omit_fields=False,
            read_header=False,
            write_header=False,
            read_records=100,
            interval=0):

    utils.expand_path(outfile)

    count = 0
    errors = 0
    duplicates = 0

    # setup the reader
    reader = rdkit_utils.create_reader(input,
                                       id_column=id_column,
                                       mol_column=mol_column,
                                       read_records=read_records,
                                       read_header=read_header,
                                       delimiter=delimiter)
    extra_field_names = reader.get_extra_field_names()

    # setup the writer
    writer = rdkit_utils.create_writer(outfile,
                                       extra_field_names=extra_field_names,
                                       delimiter=delimiter,
                                       id_column=id_column,
                                       mol_column=mol_column)

    # read the input records and write the output
    canon_smiles = {}
    id_col_type, id_col_value = utils.is_type(id_column, int)
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
                [],
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

        # deduplicate
        try:
            if mode == 'none':
                biggest = mol
                cann_smi = smi
            else:
                biggest = rdkit_utils.fragment(mol, mode)
                cann_smi = Chem.MolToSmiles(biggest)
            if cann_smi in canon_smiles:
                DmLog.emit_event("Molecule {} is duplicate of molecule {}".format(count, canon_smiles.get(cann_smi)))
                duplicates += 1
                continue
            else:
                canon_smiles[cann_smi] = count

            if omit_fields:
                for name in biggest.GetPropNames():
                    biggest.ClearProp(name)
                props = []

        except:
            errors += 1
            DmLog.emit_event('Failed to process record', count)
            # traceback.print_exc()
            continue

        # write the output
        writer.write(smi, biggest, id, props, [])

    writer.close()
    reader.close()

    return count, errors, duplicates


def main():
    # Example usage:
    #   ./rdkit_dedup.py -i data/11100.smi --outfile out.smi --delimiter tab --id-column 1 --interval 100

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='RDKit deduplicate')
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
    parser.add_argument('-m', '--mode', choices=['hac', 'mw', 'none'], default='hac',
                        help='Strategy for picking largest fragment (mw or hac or none')

    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("rdk_dedup.py: ", args)

    # special processing of delimiter to allow it to be set as a name
    delimiter = utils.read_delimiter(args.delimiter)

    t0 = time.time()
    count, errors, duplicates = process(
        args.input, args.outfile, mode=args.mode, delimiter=delimiter,
        id_column=args.id_column, mol_column=args.mol_column, omit_fields=args.omit_fields,
        read_header=args.read_header, write_header=args.write_header,
        read_records=args.read_records, interval=args.interval)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processed {} records in {} seconds. {} errors. {} duplicates'.format(
        count, duration_s, errors, duplicates))
    # Emit final 'total' cost, replacing all prior costs
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
