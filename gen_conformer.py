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

# This module is deprecated as the OpenBabel Python bindings generate some quirky structures
# that cause problems with rDock, and probably other programs too.
#
# Instead the enumerate.py module was enhanced to generate 3D structures using RDKit. Use that instead.

import os, argparse
from openbabel import pybel
import utils


def execute(infile_smiles, data_dir, interval=0):
    inputs = 0
    outputs = 0
    errors = 0
    
    with open(infile_smiles) as inf:
        for line in inf:
            inputs += 1
            
            if interval and inputs % interval == 0:
                utils.log_dm_event("Processed {} records".format(inputs))
            
            tokens1 = line.strip().split('\t')
            std_smi = tokens1[0]
            uid = tokens1[1]
            digest = tokens1[2]
            
            parts = [data_dir]
            parts.extend(utils.get_path_from_digest(digest))
            path = os.path.join(*parts)
            if not os.path.isdir(path):
                utils.log('WARNING, path', path, 'not found')
                errors += 1
                continue
                
            smi_in = os.path.join(path, digest + '.smi')
            if not os.path.exists(smi_in):
                utils.log_dm_event('WARNING, smiles file', smi_in, 'not found')
                errors += 1
                continue
                
            #utils.log('Reading', smi_in)
            with open(smi_in) as enums:
                sdf_out = os.path.join(path, digest + '.sdf')
                with pybel.Outputfile('sdf', sdf_out, overwrite=True, opt=None) as sdf:
                    count = 0
                    for l in enums:
                        outputs += 1
                        count += 1
                        tokens2 = l.strip().split('\t')
                        enum_smi = tokens2[0]
                        uid2 = tokens2[1]
                        code = tokens2[2]
                
                        mol = pybel.readstring("smi", enum_smi)
                        mol.addh()
                        mol.make3D()
                        mol.title = uid2
                        mol.data['std_smi'] = std_smi
                        mol.data['enum_smi'] = enum_smi
                        mol.data['enum_code'] = code
                        if code != 'B':
                            mol.data['parent_uuid'] = uid
                        #utils.log('Writing', std_smi, enum_smi, code)
                
                        sdf.write(mol)
            #utils.log('Wrote', digest, count)
                
    return inputs, outputs, errors



### start main execution #########################################
def main():

    # Example:
    #   python3 gen_conformer.py -i need-conf.smi --data-dir combined

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate candidates')
    parser.add_argument('-i', '--input', required=True, help="Input file as SMILES")
    parser.add_argument('--data-dir', required=True, help="Data directory")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log_dm_event("gen_conformer:", args)



    inputs, outputs, errors = execute(args.input, args.data_dir, interval=args.interval)

    utils.log_dm_event(inputs, 'inputs', outputs, 'outputs', errors, "errors")



if __name__ == "__main__":
    main()
