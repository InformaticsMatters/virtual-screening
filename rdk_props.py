#!/usr/bin/env python

# Copyright 2021 Informatics Matters Ltd.
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

import utils, mol_utils

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen


class SdfWriter:

    def __init__(self, outfile):
        self.writer = Chem.SDWriter(outfile)

    def write(self, smi, mol, id, existing_props, new_props):
        if id is not None:
            mol.SetProp('_Name', id)
        if new_props[0] is not None:
            mol.SetIntProp('hac', new_props[0])
        if new_props[1] is not None:
            mol.SetIntProp('num_rot_bonds', new_props[1])
        if new_props[2] is not None:
            mol.SetIntProp('num_rings', new_props[2])
        if new_props[3] is not None:
            mol.SetIntProp('num_aro_rings', new_props[3])
        if new_props[4] is not None:
            mol.SetIntProp('num_cc', new_props[4])
        if new_props[5] is not None:
            mol.SetIntProp('num_undef_cc', new_props[5])
        if new_props[6] is not None:
            mol.SetIntProp('num_sp3', new_props[6])

        self.writer.write(mol)

    def write_header(self, values):
        utils.log("WARNING: asked to write header for an SDF. No action will be taken.")


    def close(self):
        pass


class SmilesWriter:

    def __init__(self, outfile, sep, extra_field_names):
        self.writer = open(outfile, 'w')
        if sep is None:
            self.sep = ' '
        else:
            self.sep = sep
        self.extra_field_names = extra_field_names

    def write_header(self, values):
        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def write(self, smi, mol, id, existing_props, new_props):

        values = [smi]
        for prop in existing_props:
            if prop is not None:
                values.append(prop)
            else:
                values.append('')

        for prop in new_props:
            if prop is not None:
                values.append(str(prop))
            else:
                values.append('')


        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def close(self):
        self.writer.close()


class SdfReader:

    def __init__(self, input, id_col, recs_to_read):

        self.field_names = []
        # read a number of records to determine the field names
        if recs_to_read:
            r = self.create_reader(input)
            for i in range(recs_to_read):
                if r.atEnd():
                    break
                else:
                    mol = next(r)
                    names = mol.GetPropNames()
                    for name in names:
                        if name not in self.field_names:
                            self.field_names.append(name)

        # now create the real reader
        self.reader = self.create_reader(input)
        self.id_col = id_col

    @staticmethod
    def create_reader(input):
        if input.endswith('.gz'):
            reader = Chem.ForwardSDMolSupplier(gzip.open(input))
        else:
            reader = Chem.ForwardSDMolSupplier(input)
        return reader

    def read(self):
        try:
            props = []
            mol = next(self.reader)
            smi = Chem.MolToSmiles(mol)
            if self.id_col:
                id = mol.GetProp(self.id_col)
            else:
                id = None

            for name in self.field_names:
                if mol.HasProp(name):
                    props.append(mol.GetProp(name))
                else:
                    props.append(None)
            t = (mol, smi, id, props)
            return t

        except StopIteration:
            return None

    def get_extra_field_names(self):
        return self.field_names

    def close(self):
        pass


class SmilesReader:

    def __init__(self, input, read_header, delimiter, id_col):
        if input.endswith('.gz'):
            self.reader = gzip.open(input, 'rt')
        else:
            self.reader = open(input, 'rt')
        self.delimiter = delimiter
        if id_col is None:
            self.id_col = None
        else:
            self.id_col = int(id_col)
        self.field_names = None
        # skip header lines
        if read_header:
            line = self.reader.readline()
            tokens = self.tokenize(line)
            self.field_names = []
            for token in tokens:
                self.field_names.append(token.strip())

    def tokenize(self, line):
        line = line.strip()
        if self.delimiter is None:
            tokens = line.split()
        else:
            tokens = line.split(self.delimiter)
        stripped = []
        for token in tokens:
            stripped.append(token.strip())
        return stripped

    def read(self):
        line = self.reader.readline()
        if line:
            tokens = self.tokenize(line)
            smi = tokens[0]
            if self.id_col:
                id = tokens[self.id_col]
            else:
                id = None

            mol = Chem.MolFromSmiles(smi)
            props = []

            for i, token in enumerate(tokens):
                token = token.strip()
                if i != 0:
                    props.append(token)
                    if mol:
                        if self.field_names:
                            mol.SetProp(self.field_names[i], token)
                        else:
                            mol.SetProp('field' + str(i), token)


            t = (mol, smi, id, props)
            return t
        else:
            return None

    def get_extra_field_names(self):
        if self.field_names:
            results = []
            for i, name in enumerate(self.field_names):
                if i != 0 and i != self.id_col:
                    results.append(name)
            return results
        else:
            return None

    def close(self):
        self.reader.close()


def process(input, outfile, delimiter, id_column=None, read_header=False, write_header=False, sdf_read_records=100, interval=0):

    utils.expand_path(outfile)
    
    count = 0
    errors = 0

    # setup the reader
    if input.endswith('.sdf') or input.endswith('.sdf.gz'):
        reader = SdfReader(input, id_column, sdf_read_records)
    else:
        reader = SmilesReader(input, read_header, delimiter, id_column)
    extra_field_names = reader.get_extra_field_names()

    # setup the writer
    if outfile.endswith('.sdf') or outfile.endswith('sd'):
        writer = SdfWriter(outfile)
    else:
        writer = SmilesWriter(outfile, delimiter, extra_field_names)

    # read the input records and write the output
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break

        mol, smi, id, props = t
        if count == 0 and write_header:
            headers = ['smiles']
            if extra_field_names:
                headers.extend(extra_field_names)
            else:
                for i, prop in enumerate(props):
                    headers.append('field' + str(i + 2))
            headers.extend(['hac', 'num_rot_bonds', 'num_rings', 'num_aro_rings', 'num_cc', 'num_undef_cc',
                            'num_sp3', 'logp', 'tpsa'])
            writer.write_header(headers)

        count += 1

        if interval and count % interval == 0:
            utils.log_dm_event("Processed {} records".format(count))

        if not mol:
            errors += 1
            utils.log_dm_event("Failed to process record", count)
            continue

        # calculate the molecular props
        try:
            hac = mol.GetNumHeavyAtoms()
            num_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            num_rings = rdMolDescriptors.CalcNumRings(mol)
            num_aro_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
            num_cc, num_undef_cc = mol_utils.get_num_chiral_centers(mol)
            num_sp3 = mol_utils.get_num_sp3_centres(mol)
            logp = Crippen.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
        except:
            errors += 1
            utils.log_dm_event('Failed to process record', count)
            continue

        # write the output
        writer.write(smi, mol, id, props, (hac, num_rot_bonds, num_rings, num_aro_rings, num_cc, num_undef_cc, num_sp3, logp, tpsa))

    writer.close()
    reader.close()

    return count, errors


def main():

    # Example usage:
    #   ./rdk_props.py -i data/100000.smi -o out.sdf --interval 10000

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
    utils.log_dm_event("rdk_props.py: ", args)

    # special processing of delimiter to allow it to be set as a name

    if args.delimiter:
        if 'tab' == args.delimiter:
            delimiter = '\t'
        elif 'space' == args.delimiter:
            delimiter = None
        elif 'comma' == args.delimiter:
            delimiter = ','
        elif 'pipe' == args.delimiter:
            delimiter = '|'
        else:
            delimiter = args.delimiter
    else:
        delimiter = None

    t0 = time.time()
    count, errors = process(args.input, args.outfile, delimiter, id_column=args.id_column,
                            read_header=args.read_header, write_header=args.write_header,
                            sdf_read_records=args.sdf_read_records, interval=args.interval,)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    utils.log_dm_event('Processed {} records in {} seconds. {} errors.'.format(count, duration_s, errors))
    

if __name__ == "__main__":
    main()
