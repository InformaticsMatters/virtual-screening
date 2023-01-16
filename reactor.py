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

"""
Reaction based enumeration using RDKit.
The reaction is defined as reaction SMARTS using the --reaction argument.
The reactants are provided as SMILES files using the --reactants argument. Specify as many files as your reaction needs,
and in the appropriate order.
The products are written in SMILES format to the file specified with the --output argument. The format is:
product index reactant1 reactant2 ...
Where:
- product is the product SMILES
- index is the index (starting from 1) of the reactant combinations
- reactant1 (or reactant2 ...) is the reactant that generated that product

Note: for reactants that have multiple reactive sites there will be multiple


Example:
  python -m reactor --reaction '[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]' --reactants acids.smi amines.smi -o products.smi
"""

import argparse, time, traceback
from itertools import product

from rdkit import Chem
from rdkit.Chem import AllChem, TemplateAlign

import utils
import rdkit_utils
from dm_job_utilities.dm_log import DmLog


def read_sdf(file, index, id_field_name):
    mols = []
    count = 0
    with Chem.SDMolSupplier(file) as supplr:
        for mol in supplr:
            count += 1
            if mol:
                id = str(count)
                if id_field_name:
                    if mol.HasProp(id_field_name):
                        id = mol.GetProp(id_field_name)
                mol.SetProp('_Name', id)
                mols.append(mol)
            else:
                utils.log('Failed to read reactant', index, count)
    return mols


def read_smiles(file, index, read_header, delimiter, id_field_idx):
    mols = []
    with open(file, 'rt') as reader:
        count = 0
        if read_header:
            line = reader.readline()
        while True:
            count += 1
            line = reader.readline()
            if not line:
                break
            tokens = line.split(delimiter)
            mol = Chem.MolFromSmiles(tokens[0])
            if mol:
                id = str(count)
                if id_field_idx:
                    idx = int(id_field_idx)
                    if len(tokens) > idx:
                        id = tokens[idx]
                mol.SetProp('_Name', id)
                mols.append(mol)
            else:
                utils.log('Failed to read reactant', index, count)

    return mols


def execute(reactions, reactants, output, read_header=False, write_header=False, delimiter='\t',
            reactant_id_field=None, align=None, interval=None):

    utils.expand_path(output)

    template = None
    if align:
        template = Chem.MolFromSmiles(align)
        AllChem.Compute2DCoords(template)

    rxns = [AllChem.ReactionFromSmarts(x) for x in reactions]

    mols_arr = []
    for i, reactant in enumerate(reactants):
        DmLog.emit_event('Reading reactant', i + 1)
        if reactant.endswith('.sdf'):
            mols = read_sdf(reactant, i, reactant_id_field)
        else:
            mols = read_smiles(reactant, i, read_header, delimiter, reactant_id_field)
        mols_arr.append(mols)

    msg = 'Read reactants '
    reactant_counts = []
    for j in range(i+1):
        c = len(mols_arr[j])
        if j:
            msg += 'x'
        msg += str(c)
        reactant_counts.append(c)
    DmLog.emit_event(msg)

    num_products = 0
    num_reactions = 0
    num_combinations = 0

    field_names = ['Iteration']
    for i in range(len(reactant_counts)):
        field_names.append('Reactant' + str(i+1))
    for i in range(len(reactant_counts)):
        field_names.append('ID' + str(i+1))

    writer = rdkit_utils.create_writer(output, calc_prop_names=field_names)
    try:
        if write_header:
            headers = ['Product']
            headers.extend(field_names)
            writer.write_header(headers)

        for reacts in product(*mols_arr):
            num_combinations += 1
            reactant_smiles = [Chem.MolToSmiles(x) for x in reacts]
            reactant_ids = [x.GetProp('_Name') for x in reacts]
            product_smiles_dedup = {}
            for rxn in rxns:
                products = rxn.RunReactants(reacts)
                num_reactions += 1
                for ps in products:
                    # deduplicate the molecules using SMILES
                    for m in ps:
                        product_smiles_dedup[Chem.MolToSmiles(m)] = m

            # write the products to the output
            for smi, mol in product_smiles_dedup.items():
                props = [str(num_combinations)]
                props.extend(reactant_smiles)
                props.extend(reactant_ids)

                if template:
                    TemplateAlign.AlignMolToTemplate2D(mol, template)

                writer.write(smi, mol, str(num_combinations), [], props, smiles_prop_name='smiles')
                num_products += 1

            if interval and num_combinations % interval == 0:
                DmLog.emit_event("Processed {} reactions, generated {} products".format(num_combinations, num_products))

    finally:
        writer.close()

    return reactant_counts, num_products


def main():
    # Example:
    #   python -m reactor --reactions '[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]' --reactants data/100.smi data/100.smi -o bar.smi

    parser = argparse.ArgumentParser(description='Reactor')
    parser.add_argument('--reactions', nargs='+', required=True, help="Reaction SMARTS")
    parser.add_argument('--reactants', nargs='+', required=True, help="Reactant files as .smi")
    parser.add_argument('-o', '--output', required=True, help="Output file (.smi or .sdf)")
    parser.add_argument('--read-header', action='store_true', help="Reactant files have header line")
    parser.add_argument('--write-header', action='store_true', help="Add header line to output")
    parser.add_argument('--delimiter', help="Delimiter when using SMILES")
    parser.add_argument('--reactant-id-field', help="Index or name of the IDs of the reactants")
    parser.add_argument('--align', help="Align molecules to this core structure (SMILES). SDF output only.")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("reactor: ", args)

    delimiter = utils.read_delimiter(args.delimiter)

    start = time.time()
    reactant_counts, product_count = execute(args.reactions, args.reactants, args.output,
                                             read_header=args.read_header, write_header=args.write_header,
                                             delimiter=delimiter, reactant_id_field=args.reactant_id_field,
                                             align=args.align, interval=args.interval)
    end = time.time()

    msg = 'Reactants: '
    q = 1
    for i, c in enumerate(reactant_counts):
        q = q * c
        if i:
            msg += 'x'
        msg += str(c)
    DmLog.emit_event(msg, 'Reactions:', q, 'Products:', str(product_count),  'Time (s):', end - start)
    DmLog.emit_cost(product_count)


if __name__ == "__main__":
    main()
