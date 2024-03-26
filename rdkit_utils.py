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

import csv
import gzip
from rdkit import Chem
import utils

ID_COL_NAME = "ID"
SMILES_COL_NAME = "SMILES"


class SdfWriter:

    def __init__(self, outfile, prop_names):
        if outfile.endswith('.gz'):
            self.gzip = gzip.open(outfile, 'wt')
            self.writer = Chem.SDWriter(self.gzip)
        else:
            self.gzip = None
            self.writer = Chem.SDWriter(outfile)
        self.prop_names = prop_names

    def write(self, smi, mol, id, existing_props, new_props, smiles_prop_name=None):
        if not mol:
            mol = Chem.MolFromSmiles(smi)
        if id is not None:
            mol.SetProp('_Name', id)
        if smiles_prop_name is not None:
            mol.SetProp(smiles_prop_name, smi)

        for i, prop_name in enumerate(self.prop_names):
            value = new_props[i]
            if prop_name is not None:
                mol.SetProp(prop_name, str(value))

        self.writer.write(mol)

    def write_header(self, values):
        utils.log("INFO: asked to write header for an SDF. No action will be taken.")

    def close(self):
        self.writer.close()
        if self.gzip:
            self.gzip.close()


class SmilesWriter:

    def __init__(self, outfile, sep, extra_field_names, id_column=None, mol_column=0):
        self.writer = open(outfile, 'w')
        if sep is None:
            self.sep = ' '
        else:
            self.sep = sep
        self.extra_field_names = extra_field_names


        self.id_column = None if id_column is None else int(id_column)
        self.mol_column = None if mol_column is None else int(mol_column)

    def write_header(self, values):
        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def write(self, smi, mol, id, existing_props, new_props, smiles_prop_name=None):
        if type(id) == str and self.sep in id:
            id = '"' + id + '"'

        if self.id_column is not None:
            if self.id_column < self.mol_column:
                values = [id, smi]
            else:
                values = [smi, id]
        else:
            values = [smi]

        for prop in existing_props:
            if prop is not None:
                if type(prop) == str and self.sep in prop:
                    values.append('"' + prop + '"')
                else:
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
        if id_col == "_Name":
            self.field_names.append(ID_COL_NAME)
        # read a number of records to determine the field names
        if recs_to_read:
            r = self.create_reader(input)
            for i in range(recs_to_read):
                if r.atEnd():
                    break
                else:
                    try:
                        mol = next(r)
                        if mol:
                            names = mol.GetPropNames()
                            for name in names:
                                if name not in self.field_names:
                                    self.field_names.append(name)
                    except StopIteration:
                        break

        # now create the real reader
        self.reader = self.create_reader(input)
        self.id_col = id_col

    def get_mol_field_name(self):
        return None

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
                if mol.HasProp(self.id_col):
                    id = mol.GetProp(self.id_col)
                else:
                    id = None
            else:
                id = None
            if self.id_col == "_Name":
                props.append(id)

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

    def __init__(self, input, read_header, delimiter, id_column, mol_column, recs_to_read):
        self.delimiter = delimiter if delimiter is not None else ' '
        if mol_column is None:
            if id_column == 0:
                self.mol_column = 1
            else:
                self.mol_column = 0
        else:
            self.mol_column = mol_column

        if id_column is None:
            self.id_column = None
        else:
            self.id_column = int(id_column)

        tmp_reader, file_reader = self.create_readers(input)
        # read header line
        if read_header:
            tokens = next(tmp_reader)
            # tokens = self.tokenize(line)
            self.field_names = []
            for token in tokens:
                self.field_names.append(token.strip())
        else:
            self.field_names = []
            for i in range(0, max(1, self.mol_column + 1, 0 if self.id_column is None else self.id_column + 1)):
                self.field_names.append(None)
            self.field_names[self.mol_column] = SMILES_COL_NAME
            if self.id_column is not None:
                self.field_names[self.id_column] = ID_COL_NAME

        max_num_tokens = 0
        for i in range(0, recs_to_read):
            line = next(tmp_reader, None)
            if not line:
                break
            else:
                num_tokens = len(line)
                if num_tokens > max_num_tokens:
                    max_num_tokens = num_tokens

        if max_num_tokens > len(self.field_names):
            for i in range(len(self.field_names), max_num_tokens):
                self.field_names.append('field' + str(i + 1))

        file_reader.close()

        # now create the read reader and discard the header
        self.reader, self.file = self.create_readers(input)
        if read_header:
            line = next(self.reader)

    def create_readers(self, input):
        if input.endswith('.gz'):
            r = gzip.open(input, 'rt')
            return csv.reader(r, delimiter=self.delimiter), r
        else:
            r = open(input, 'rt')
            return csv.reader(r, delimiter=self.delimiter), r

    def get_mol_field_name(self):
        if self.field_names:
            return self.field_names[self.mol_column]
        else:
            return None

    def read(self):
        tokens = next(self.reader, None)
        if tokens:
            smi = tokens[self.mol_column]
            if self.id_column is not None:
                id = tokens[self.id_column]
            else:
                id = None

            mol = Chem.MolFromSmiles(smi)
            props = []

            for i, token in enumerate(tokens):
                token = token.strip()
                if not (i == self.mol_column or i == self.id_column):
                    props.append(token)
                    if mol:
                        if self.field_names and len(self.field_names) > i:
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
                if i != 0 and i != self.id_column:
                    results.append(name)
            return results
        else:
            return []

    def close(self):
        self.file.close()


def generate_headers(
        id_col_type,
        id_col_value,
        mol_field_name: str,
        field_names: list[str],
        calc_prop_names:
        list[str],
        omit_fields: bool):
    """
    Generate the headers for when writing a tab or comma separated files
    
    :param id_col_type: The type of ID column that was specified.
        -1 String for the SDF field name
        +1 int for column index for TAB
        0 for None (no ID col)
    :param id_col_value: The value specified for the ID column
        -1: the string specified
        +1: the value specified as an int
        0: None
    :param mol_field_name: The name of the mol column, or None for SDF
    :param field_names: The names of the fields that were in the input that will be added to the input
    :param calc_prop_names: The names of the new fields to be added
    :param omit_fields: Do not add the input fields (except for the mol and ID)
    :return: List of the header names 
    """
    headers = []
    mol_header = mol_field_name if mol_field_name else SMILES_COL_NAME

    if id_col_type == 0:  # id_col was None
        if omit_fields or mol_header not in field_names:
            headers.append(mol_header)
    elif id_col_type == 1:  # was an int so CSV
        if omit_fields or mol_header not in field_names:
            headers.append(mol_header)
        if omit_fields:
            headers.append(field_names[id_col_value])
    else:  # id_col was string, so SDF field name or _Name
        headers.append(mol_field_name if mol_field_name else SMILES_COL_NAME)
        if id_col_value == '_Name':
            id_name = ID_COL_NAME
        else:
            id_name = id_col_value
        if id_name not in field_names or omit_fields:
            headers.append(id_name)

    if not omit_fields:
        headers.extend(field_names)

    headers.extend(calc_prop_names)
    return headers


def create_reader(input, type=None, id_column=None, mol_column=None, read_records=100, read_header=False, delimiter='\t'):
    if type is None:
        if input.endswith('.sdf') or input.endswith('.sdf.gz') or input.endswith('.sd') or input.endswith('.sd.gz'):
            type = 'sdf'
        else:
            type = 'smi'

    if type == 'sdf':
        return SdfReader(input, id_column, read_records)
    elif type == 'smi':
        return SmilesReader(input, read_header, delimiter, id_column, mol_column, read_records)
    else:
        raise ValueError('Unexpected file type', type)


def create_writer(outfile, delimiter='\t', extra_field_names=[], calc_prop_names=[], id_column=None, mol_column=0):
    if outfile.endswith('.sdf') or outfile.endswith('sd'):
        return SdfWriter(outfile, calc_prop_names)
    else:
        return SmilesWriter(outfile, delimiter, extra_field_names, id_column=id_column, mol_column=mol_column)


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
                    utils.log("ERROR! " + str(lines[0]) + "unknown charge flag: " + str(chg[1]))
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


def rdk_read_mols(input):
    """
    Read molecules from a single file (.mol or .sdf)
    """
    if input.endswith('.mol'):
        mol = Chem.MolFromMolFile(input)
        return [mol]
    else:
        supplr = rdk_mol_supplier(input)
        mols = [m for m in supplr]
        return mols


def rdk_read_molecule_files(inputs):
    """
    Read input molecules.
    A list of inputs is specified, each one containing either a single filename or a comma separated list of filenames
    (no spaces). The files can either be .mol files with a single molecule or .sdf files with multiple molecules.
    The molecules in all the files specified are read and returned as a list of molecules.
    Examples:
    ['mols.sdf']                          - single SDF with one or more molecules
    ['mol1.mol', 'mol2.mol', 'mol3.mol']  - 3 molfiles as separate elements of the list
    ['mol1.mol,mol2.mol,mol3.mol']        - 3 molfiles as a single comma separated element of the list
    ['mol1.mol,mol2.mol', 'mols.sdf']     - 3 molfiles plus 1 SDF

    :param inputs: Input filenames
    :return: List of molecules that have been read
    """
    mols = []
    for input in inputs:
        tokens = input.split(',')
        for token in tokens:
            if token.endswith('.mol'):
                m = Chem.MolFromMolFile(token)
                if m:
                    mols.append(m)
                else:
                    DmLog.emit_event('WARNING: could not process', token)
            else:
                ms = rdk_read_mols(token)
                for i, m in enumerate(ms):
                    if m:
                        mols.append(m)
                    else:
                        DmLog.emit_event('WARNING: could not process molecule', i, 'from', token)
    return mols


def rdk_merge_mols(inputs):
    """
    Merge multiple molecules into a single molecule
    :param input:
    :return:
    """
    merged_mol = Chem.RWMol()
    mols = rdk_read_molecule_files(inputs)
    for i, mol in enumerate(mols):
        if mol:
            merged_mol.InsertMol(mol)
        else:
            DmLog.emit_event('WARNING: could not process molecule', i)
    Chem.SanitizeMol(merged_mol)
    return merged_mol, i + 1


def rdk_mol_supplier(input):
    """
    Generate a ForwardSDMolSupplier from a .sdf or .sdf.gz file
    """
    if input.endswith('.sdf'):
        suppl = Chem.ForwardSDMolSupplier(input)
    elif input.endswith('.sdf.gz'):
        gz = gzip.open(input, 'rb')
        suppl = Chem.ForwardSDMolSupplier(gz)
    else:
        raise ValueError('Unsupported file type. Must be .sdf or .sdf.gz. Found ' + input)
    return suppl


def fragment(mol, mode):
    """
    Generate the largest fragment in the molecule e.g. typically a desalt operation
    :param mol: The molecule to fragment
    :param mode: The strategy for picking the largest (mw or hac)
    :return:
    """
    frags = Chem.GetMolFrags(mol, asMols=True)

    if len(frags) == 1:
        return mol
    else:
        # TODO - handle ties
        biggest_index = -1
        i = 0
        if mode == 'hac':
            biggest_count = 0
            for frag in frags:
                hac = frag.GetNumHeavyAtoms()
                if hac > biggest_count:
                    biggest_count = hac
                    biggest_mol = frag
                    biggest_index = i
                i+=1
            # utils.log("Chose fragment", biggest_index, "from", len(frags), "based on HAC")
        elif mode == 'mw':
            biggest_mw = 0
            for frag in frags:
                mw = Descriptors.MolWt(frag)
                if mw > biggest_mw:
                    biggest_mw = mw
                    biggest_mol = frag
                    biggest_index = i
                i+=1
            # utils.log("Chose fragment", biggest_index, "from", len(frags), "based on MW")
        else:
            raise ValueError('Invalid fragment mode:',mode)

        # copy the properties across
        for name in mol.GetPropNames():
            biggest_mol.SetProp(name, mol.GetProp(name))

        # _Name is a magical property that is not in the ones returned by GetPropNames
        if '_Name' in mol.GetPropNames():
            biggest_mol.SetProp("_Name", mol.GetProp("_Name"))

        return biggest_mol


def fragmentAndFingerprint(reader, mols, data, fps, descriptor, fragmentMethod='hac', outputFragment=False):
    """
    Fragment the molecule if it has multiple fragments and generate fingerprints on the fragment.

    :param reader: MolSupplier from which to read the molecules
    :param mols: List to which the molecules are added
    :param fps: List to which the fingerprints are added
    :param descriptor: Function to generate the fingerprints from the molecule
    :param fragmentMethod: The fragmentation method to use when there are multiple fragments (hac or mw)
    :param outputFragment: Boolean that specifies whether to write the fragment or the original molecule to the mols list
    :return: The number of errors encountered
    """
    errors = 0
    count = 0
    while True:
        count += 1
        t = reader.read()
        # break if no more data to read
        if not t:
            break
        mol, smi, id, props = t
        if not mol:
            utils.log('Failed to read molecule', count)
            errors += 1
            continue

        frag = fragment(mol, fragmentMethod)
        d = descriptor(frag)
        if d:
            if outputFragment:
                mols.append(frag)
                data.append((id, Chem.MolToSmiles(frag), props))
            else:
                mols.append(mol)
                data.append((id, Chem.MolToSmiles(mol), props))
            fps.append(d)
            continue
    return errors


def check_molecules_are_3d(input, num_to_check=10):
    if input.endswith('.gz'):
        suppl = Chem.ForwardSDMolSupplier(gzip.open(input))
    else:
        suppl = Chem.ForwardSDMolSupplier(input)

    count = 0
    for mol in suppl:
        conf = mol.GetConformer()
        if conf is None or not conf.Is3D():
            return False
        count += 1
        if count == num_to_check:
            break
    return True
