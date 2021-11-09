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
from rdkit.Chem import rdMolDescriptors

class SdfWriter():

    def __init__(self, outfile):
        self.writer = Chem.SDWriter(outfile)

    def write(self, smi, mol, id, props):
        if id is not None:
            mol.SetProp('_Name', id)
        if props[0] is not None:
            mol.SetIntProp('hac', props[0])
        if props[1] is not None:
            mol.SetIntProp('num_rot_bonds', props[1])
        if props[2] is not None:
            mol.SetIntProp('num_rings', props[2])
        if props[3] is not None:
            mol.SetIntProp('num_aro_rings', props[3])
        if props[4] is not None:
            mol.SetIntProp('num_cc', props[4])
        if props[5] is not None:
            mol.SetIntProp('num_undef_cc', props[5])
        if props[6] is not None:
            mol.SetIntProp('num_sp3', props[6])

        self.writer.write(mol)

    def close(self):
        pass

class SmilesWriter:

    def __init__(self, outfile, sep, extra_field_names):
        self.writer = open(outfile, 'w')
        self.sep = sep
        self.extra_field_names = extra_field_names

    def write_header(self, values):
        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def write(self, smi, mol, id, props):

        values = [smi, id]
        if self.extra_field_names:
            for name in self.extra_field_names:
                if mol.HasProp(name):
                    values.append(mol.GetProp(name))
                else:
                    values.append('')

        for prop in props:
            if prop is None:
                values.append('')
            else:
                values.append(str(prop))

        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def close(self):
        self.writer.close()


class SdfReader:

    def __init__(self, input, id_col, recs_to_read):

        self.field_names = set()
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
                        self.field_names.add(name)

        # now create the real reader
        self.reader = self.create_reader(input)
        self.id_col = id_col

    @staticmethod
    def create_reader(input):
        if input.endswith('.gz'):
            reader = Chem.ForwardSDMolSupplier(gzip.open(input))
        else:
            reader = Chem.SDMolSupplier(input)
        return reader

    def read(self):
        if self.reader.atEnd():
            return None
        else:
            mol = next(self.reader)
            smi = Chem.MolToSmiles(mol)
            if self.id_col:
                id = mol.GetProp(self.id_col)
            else:
                id = None
            t = (mol, smi, id)
            return t

    def get_extra_field_names(self):
        return self.field_names

    def close(self):
        pass


class SmilesReader():

    def __init__(self, input, skip_lines, delimiter, id_col):
        if input.endswith('.gz'):
            self.reader = gzip.open(input, 'rt')
        else:
            self.reader = open(input, 'rt')
        self.delimiter = delimiter
        self.id_col = int(id_col)
        self.field_names = None
        # skip header lines
        if skip_lines:
            for i in range(skip_lines):
                line = self.reader.readline()
                if i == 0:
                    tokens = self.tokenize(line)
                    self.field_names = []
                    for token in tokens:
                        self.field_names.append(token.strip())

    def tokenize(self, line):
        line = line.strip()
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
            if mol:
                for i, token in enumerate(tokens):
                    if i != 0 and i != self.id_col and self.field_names:
                        mol.SetProp(self.field_names[i], token)

            t = (mol, smi, id)
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


def process(input, outfile, delimiter, name_column=None, no_header=False, skip_lines=0, sdf_read_records=100, interval=0):

    utils.expand_path(outfile)
    
    count = 0
    errors = 0

    # setup the reader
    if input.endswith('.sdf') or input.endswith('.sdf.gz'):
        reader = SdfReader(input, name_column, sdf_read_records)
    else:
        reader = SmilesReader(input, skip_lines, delimiter, name_column)
    extra_field_names = reader.get_extra_field_names()

    # setup the writer
    if outfile.endswith('.smi') or outfile.endswith('txt'):
        writer = SmilesWriter(outfile, "\t", extra_field_names)
        if not no_header:
            headers = ['smiles', 'id']
            if extra_field_names:
                headers.extend(extra_field_names)
            headers.extend(['hac', 'num_rot_bonds', 'num_rings', 'num_aro_rings', 'num_cc', 'num_undef_cc', 'num_sp3'])
            writer.write_header(headers)
    elif outfile.endswith('.sdf'):
        writer = SdfWriter(outfile)
    else:
        raise ValueError('Outfile must be .sdf .smi or .txt')

    # read the input records and write the output
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break

        count += 1
        mol, smi, id = t

        if interval and count % interval == 0:
            utils.log_dm_event("Processed {} records".format(count))

        if not mol:
            errors += 1
            utils.log_dm_event("Failed to process record", count)
            continue

        if not id:
            id = str(count)

        # calculate the molecular props
        try:
            hac = mol.GetNumHeavyAtoms()
            num_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            num_rings = rdMolDescriptors.CalcNumRings(mol)
            num_aro_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
            num_cc, num_undef_cc = mol_utils.get_num_chiral_centers(mol)
            num_sp3 = mol_utils.get_num_sp3_centres(mol)
        except:
            errors += 1
            utils.log_dm_event('Failed to process record', count)
            continue

        # write the output
        writer.write(smi, mol, id, (hac, num_rot_bonds, num_rings, num_aro_rings, num_cc, num_undef_cc, num_sp3))

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
    parser.add_argument('-d', '--delimiter', default='\t', help="Delimiter when using SMILES")
    parser.add_argument('-n', '--name-column', help="Column for name field (zero based integer for .smi, text for SDF)")
    parser.add_argument('--skip-lines', default=0, type=int, help="Skip this many lines e.g. use 1 for skipping a header line")
    parser.add_argument('--no-header', action='store_true', help='Do not write header line when writing .smi or .txt')
    parser.add_argument('--sdf-read-records', default=100, type=int, help="Read this many SDF records to determine field names")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log_dm_event("rdk_props.py: ", args)
    
    t0 = time.time()
    count, errors = process(args.input, args.outfile, args.delimiter, name_column=args.name_column, interval=args.interval,
                            skip_lines=args.skip_lines, no_header=args.no_header, sdf_read_records=args.sdf_read_records)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    utils.log_dm_event('Processed {} records in {} seconds. {} errors.'.format(count, duration_s, errors))
    

if __name__ == "__main__":
    main()
