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

import utils, rdkit_utils
from dm_job_utilities.dm_log import DmLog

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen


def process(input, outfile, delimiter, id_column=None, read_header=False, write_header=False, sdf_read_records=100, interval=0):

    utils.expand_path(outfile)
    
    count = 0
    errors = 0

    calc_prop_names = ['hac', 'num_rot_bonds', 'num_rings', 'num_aro_rings', 'num_cc', 'num_undef_cc',
                       'num_sp3', 'logp', 'tpsa']

    # setup the reader
    reader = rdkit_utils.create_reader(input, id_column=id_column, sdf_read_records=sdf_read_records,
                                       read_header=read_header, delimiter=delimiter)
    extra_field_names = reader.get_extra_field_names()

    # setup the writer
    writer = rdkit_utils.create_writer(outfile, extra_field_names=extra_field_names, calc_prop_names=calc_prop_names,
                                       delimiter=delimiter)

    # read the input records and write the output
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break

        mol, smi, id, props = t
        if count == 0 and write_header:
            headers = rdkit_utils.generate_header_values(extra_field_names, calc_prop_names)
            writer.write_header(headers)

        count += 1

        if interval and count % interval == 0:
            DmLog.emit_event("Processed {} records".format(count))
        if count % 10000 == 0:
            # Emit a 'total' cost, replacing all prior costs
            DmLog.emit_cost(count)

        if not mol:
            errors += 1
            DmLog.emit_event("Failed to process record", count)
            continue

        # calculate the molecular props
        try:
            hac = mol.GetNumHeavyAtoms()
            num_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            num_rings = rdMolDescriptors.CalcNumRings(mol)
            num_aro_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
            num_cc, num_undef_cc = rdkit_utils.get_num_chiral_centers(mol)
            num_sp3 = rdkit_utils.get_num_sp3_centres(mol)
            logp = Crippen.MolLogP(mol)
            # logp values have a silly number of decimal places
            if logp is not None:
                logp = "%.2f" % logp
            tpsa = rdMolDescriptors.CalcTPSA(mol)
        except:
            errors += 1
            DmLog.emit_event('Failed to process record', count)
            continue

        # write the output
        writer.write(smi, mol, id, props, (hac, num_rot_bonds, num_rings, num_aro_rings, num_cc, num_undef_cc, num_sp3, logp, tpsa))

    writer.close()
    reader.close()

    return count, errors


def main():

    # Example usage:
    #   ./rdk_props.py -i data/100000.smi --outfile out.sdf --delimiter tab --interval 10000

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Calc RDKit props')
    parser.add_argument('-i', '--input', required=True, help="Input file as SMILES or SDF")
    parser.add_argument('-o', '--outfile', required=True, help="Output file as SMILES or SDF")
    # to pass tab as the delimiter specify it as $'\t' or use one of the symbolic names 'comma', 'tab', 'space' or 'pipe'
    parser.add_argument('-d', '--delimiter', help="Delimiter when using SMILES")
    parser.add_argument('--id-column', help="Column for name field (zero based integer for .smi, text for SDF)")
    parser.add_argument('--read-header', action='store_true', help="Read a header line with the field names when reading .smi or .txt")
    parser.add_argument('--write-header', action='store_true', help='Write a header line when writing .smi or .txt')
    parser.add_argument('--sdf-read-records', default=100, type=int, help="Read this many SDF records to determine field names")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("rdk_props.py: ", args)

    # special processing of delimiter to allow it to be set as a name

    delimiter = utils.read_delimiter(args.delimiter)

    t0 = time.time()
    count, errors = process(args.input, args.outfile, delimiter, id_column=args.id_column,
                            read_header=args.read_header, write_header=args.write_header,
                            sdf_read_records=args.sdf_read_records, interval=args.interval,)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processed {} records in {} seconds. {} errors.'.format(count, duration_s, errors))
    # Emit final 'total' cost, replacing all prior costs
    DmLog.emit_cost(count)
    

if __name__ == "__main__":
    main()
