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

import gzip
from rdkit import Chem
import utils

class SdfWriter:

    def __init__(self, outfile, prop_names):
        self.writer = Chem.SDWriter(outfile)
        self.prop_names = prop_names

    def write(self, smi, mol, id, existing_props, new_props):
        if id is not None:
            mol.SetProp('_Name', id)

        for i, prop_name in enumerate(self.prop_names):
            value = new_props[i]
            if prop_name is not None:
                mol.SetProp(prop_name, str(value))

        self.writer.write(mol)

    @staticmethod
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


def create_reader(input, type=None, id_column=None, sdf_read_records=100, read_header=False, delimiter='\t'):
    if type is None:
        if input.endswith('.sdf') or input.endswith('.sdf.gz') or input.endswith('.sd') or input.endswith('.sd.gz'):
            type = 'sdf'
        else:
            type = 'smi'

    if type == 'sdf':
        return SdfReader(input, id_column, sdf_read_records)
    elif type == 'smi':
        return SmilesReader(input, read_header, delimiter, id_column)
    else:
        raise ValueError('Unexpected file type', type)


def updateChargeFlagInAtomBlock(mb):
    """
    Add data for the charges to the atom block. This data is now deprecated and should be specified using "M  CHG" lines
    but some old software such as rDock only handle the old syntax. RDKit only supports the new syntax so this method
    handles adding the old syntax as well.
    This function is based on work by Jose Manuel Gally that can be found here:
    See https://sourceforge.net/p/rdkit/mailman/message/36425493/
    """
    f="{:>10s}"*3+"{:>2}{:>4s}"+"{:>3s}"*11
    chgs = []    # list of charges
    lines = mb.split("\n")
    #if mb[0] == '' or mb[0] == "\n":
    #    del lines[0]
    CTAB = lines[3]
    atomCount = int(CTAB.split()[0])
    # parse mb line per line
    for l in lines:
        # look for M CHG property
        if l[0:6] == "M  CHG":
            records = l.split()[3:]    # M  CHG X is not needed for parsing, the info we want comes afterwards
            # record each charge into a list
            for i in range(0,len(records),2):
                idx = records[i]
                chg = records[i+1]
                chgs.append((int(idx), int(chg)))    # sort tuples by first element?
            break    # stop iterating

    # sort by idx in order to parse the molblock only once more
    chgs = sorted(chgs, key=lambda x: x[0])

    # that we have a list for the current molblock, attribute each charges
    for chg in chgs:
        i=4
        while i < 4+atomCount:    # do not read from beginning each time, rather continue parsing mb!
            # when finding the idx of the atom we want to update, extract all fields and rewrite whole sequence
            if i-3 == chg[0]:    # -4 to take into account the CTAB headers, +1 because idx begin at 1 and not 0
                fields = lines[i].split()
                x=fields[0]
                y=fields[1]
                z=fields[2]
                symb=fields[3]
                massDiff=fields[4]
                charge=fields[5]
                sp=fields[6]
                hc=fields[7]
                scb=fields[8]
                v=fields[9]
                hd=fields[10]
                nu1=fields[11]
                nu2=fields[12]
                aamn=fields[13]
                irf=fields[14]
                ecf=fields[15]
                # update charge flag
                if chg[1] == -1:
                    charge = '5'
                elif chg[1] == -2:
                    charge = '6'
                elif chg[1] == -3:
                    charge = '7'
                elif chg[1] == 1:
                    charge = '3'
                elif chg[1] == 2:
                    charge = '2'
                elif chg[1] == 3:
                    charge = '1'
                else:
                    print("ERROR! " + str(lines[0]) + "unknown charge flag: " + str(chg[1]))    # print name then go to next chg
                    break
                # update modatom block line
                lines[i] = f.format(x,y,z,symb,massDiff,charge,sp,hc,scb,v,hd,nu1,nu2,aamn,irf,ecf)
            i+=1
    upmb = "\n".join(lines)
    return(upmb)


def get_num_chiral_centers(mol):
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    undef_cc = 0
    for cc in chiral_centers:
        if cc[1] == '?':
            undef_cc += 1
    return len(chiral_centers), undef_cc


def get_num_sp3_centres(mol):
    return sum((x.GetHybridization() == Chem.HybridizationType.SP3) for x in mol.GetAtoms())


def sdf_read_mol(input):
    txt = sdf_read_block(input)
    if txt == None:
        return None
    mol = Chem.MolFromMolBlock(txt)
    print('NumProps:', mol.GetNumProps())
    return (txt, mol)


def sdf_record_gen(hnd):
    """A generator for text records fom a SD file
    """
    mol_text_tmp = ""
    while 1:
        line = hnd.readline()
        if not line:
            return
        line = line.decode("utf-8")
        mol_text_tmp += line
        if line.startswith("$$$$"):
            mol_text = mol_text_tmp
            mol_text_tmp = ""
            yield mol_text


def rdk_read_single_mol(input):
    # read the molecule. Can be SDF or Mol format
    # if SDF then the first molecule is used.
    if input.endswith('.mol'):
        mol = Chem.MolFromMolFile(input)
    elif input.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(input)
        mol = next(suppl)
    elif input.endswith('.sdf.gz'):
        with gzip.open(input, 'rb') as gz:
            suppl = Chem.ForwardSDMolSupplier(gz)
            mol = next(suppl)
    else:
        raise ValueError('Unsupported file type. Must be .mol .sdf or .sdf.gz. Found ' + input)
    return mol


def rdk_merge_mols(inputs):
    """
    Merge multiple molecules into a single molecule
    :param input:
    :return:
    """
    merged_mol = Chem.RWMol()
    count = 0
    for input in inputs:
        if input.endswith('.mol'):
            mol = Chem.MolFromMolFile(input)
            merged_mol.InsertMol(mol)
            count += 1
        elif input.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(input)
            for mol in suppl:
                merged_mol.InsertMol(mol)
                count += 1
        elif input.endswith('.sdf.gz'):
            with gzip.open(input, 'rb') as gz:
                suppl = Chem.ForwardSDMolSupplier(gz)
                for mol in suppl:
                    merged_mol.InsertMol(mol)
                    count += 1
        else:
            raise ValueError('Unsupported file type. Must be .mol .sdf or .sdf.gz. Found ' + input)
    Chem.SanitizeMol(merged_mol)
    return merged_mol, count


def rdk_mol_supplier(input):
    if input.endswith('.sdf'):
        suppl = Chem.ForwardSDMolSupplier(input)
    elif input.endswith('.sdf.gz'):
        gz = gzip.open(input, 'rb')
        suppl = Chem.ForwardSDMolSupplier(gz)
    else:
        raise ValueError('Unsupported file type. Must be .mol .sdf or .sdf.gz. Found ' + input)
    return suppl
