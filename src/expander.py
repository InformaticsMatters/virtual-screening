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

# Cypher query that is executed:
# MATCH p=(m:F2)-[:FRAG*1..2]-(e:Mol)<-[:NonIso*0..1]-(c:Mol)
# WHERE m.smiles=$smiles AND e.smiles <> $smiles
# RETURN p LIMIT $limit
#
"""
Expand a set of molecules using the fragment network
"""

from rdkit import Chem
from rdkit import RDLogger

import os, argparse, json
import requests
from collections import OrderedDict

import  utils
from standardize_molecule import standardize_to_noniso_smiles

RDLogger.logger().setLevel(RDLogger.ERROR)

def process(input, outdir='./', hac_min=None, hac_max=None, rac_min=None, rac_max=None, hops=1,
            server='https://fragnet-search.xchem-dev.diamond.ac.uk/fragnet-search-api',
            token=None, index_as_filename=False, exclude_inputs=False, excludes=None):

    url = server + '/fragnet-search/rest/v2/search/expand/'
    params = {'hops': hops, 'hac_min': hac_min, 'hac_max': hac_max, 'rac_min': rac_min, 'rac_max': rac_max}
    # utils.log('URL', url, "Params:" + str(params))
    headers = {'Content-Type': 'chemical/x-daylight-smiles'}
    if token:
        headers['Authorization'] = 'bearer ' + token

    molecules = []
    duplicates = set()
    excluded_mols = set()

    # process the input, keeping only non-duplicates
    count = 0
    num_duplicates = 0
    supplr = Chem.SDMolSupplier(input)
    for mol in supplr:
        count += 1
        name = mol.GetProp('_Name')
        smiles = mol.GetProp('std_smi')
        std_smi, mol = standardize_to_noniso_smiles(smiles)
        if std_smi in duplicates:
            num_duplicates += 1
            utils.log('Duplicate:', std_smi)
            continue
        else:
            molecules.append((count, std_smi, mol, name))
            if exclude_inputs:
                excluded_mols.add(std_smi)

    # process the molecules to exclude
    if excludes:
        for excl in excludes:
            supplr = Chem.SDMolSupplier(excl)
            for mol in supplr:
                smiles = mol.GetProp('std_smi')
                excluded_mols.add(smiles)

    for item in molecules:
        count = item[0]
        std_smi = item[1]
        mol = item[2]
        name = item[3]
        utils.log('Processing', count, name, std_smi)

        r = requests.get(url + requests.utils.quote(std_smi), params=params, headers=headers)
        if r.status_code == requests.codes.ok:
            j = r.json()
            num_mols = j['size']
            utils.log('JSON OK', num_mols, 'results')
            basename = None
            if not index_as_filename:
                if mol.HasProp('_Name'):
                    name = mol.GetProp('_Name')
                    if name:
                        basename = name
            if not basename:
                basename = str(count)
            
            if not outdir.endswith('/'):
                outdir = outdir + '/'    
            utils.expand_path(outdir)

            with open(outdir + basename + '.json', 'w') as f:
                s = json.dumps(j)
                f.write(s)

            with open(outdir + basename + '.smi', 'w') as f:
                members = j['members']
                incl = 0
                excl = 0
                for member in members:
                    if member['smiles'] in excluded_mols:
                        excl += 1
                    else:
                        incl += 1
                        f.write(member['smiles'] + '\n')
                utils.log('Included', incl, 'excluded', excl)
        else:
            utils.log('Request failed with status code', r.status_code, r.text)
    return count, num_duplicates


def main():
    # Example usage:
    # 1. Create keycloak token:
    # export KEYCLOAK_TOKEN=$(curl -d "grant_type=password" -d "client_id=fragnet-search-ui" -d "username=<username>" -d "password=<password>" \
    #   https://squonk.it/auth/realms/squonk/protocol/openid-connect/token 2> /dev/null | jq -r '.access_token')
    #
    # 2. Run the module:
    #  python expander.py -i inputs.smi --hac-min 3 --hac-max 3 --rac-min 1 --rac-max 1 --hops 1 --token $KEYCLOAK_TOKEN

    parser = argparse.ArgumentParser(description='Fragnet expand')
    parser.add_argument('-i', '--input', help='Input SDF file')
    parser.add_argument('-o', '--outdir', default='./', help='Directory for outputs')
    parser.add_argument('--hac-min', type=int, default=3, help='The min change in heavy atom count')
    parser.add_argument('--hac-max', type=int, default=3, help='The max change in heavy atom count')
    parser.add_argument('--rac-min', type=int, default=1, help='The min change in ring atom count')
    parser.add_argument('--rac-max', type=int, default=1, help='The max change in ring atom count')
    parser.add_argument('--hops', type=int, default=1, help='The number of graph traversals (hops)')
    parser.add_argument('-s', '--server', default='https://fragnet-search.xchem-dev.diamond.ac.uk/fragnet-search-api', help='The fragnet search server')
    parser.add_argument('--token', help='Keycloak auth token (or specify as KEYCLOAK_TOKEN env variable')
    parser.add_argument('--index-as-filename', action='store_true', help='Use the molecule index as the file name instead of the molecule name')
    parser.add_argument('--exclude-inputs', action='store_true', help='Exclude inputs in the .smi output')
    parser.add_argument('-e', '--excludes', nargs='*', help='Exclude these from the output')

    args = parser.parse_args()
    utils.log("FragnetExpand Args: ", args)

    if args.token:
        auth_token = args.token
    else:
        auth_token = os.getenv('KEYCLOAK_TOKEN')
    if not auth_token:
        utils.log('WARNING: no authentication token found in environment variable KEYCLOAK_TOKEN')

    # this does the processing
    count, duplicates = process(args.input, outdir=args.outdir, hac_min=args.hac_min, hac_max=args.hac_max, rac_min=args.rac_min, rac_max=args.rac_max,
            hops=args.hops,
            server=args.server, token=auth_token,
            index_as_filename=args.index_as_filename,
            exclude_inputs=args.exclude_inputs,
            excludes=args.excludes)

    utils.log('Processed', count, 'molecules,', duplicates, 'duplicates')


if __name__ == "__main__":
    main()
