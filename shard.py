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


import argparse, os, json, gzip, time, uuid
import utils
from hashlib import sha256
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from standardize_molecule import standardize_to_iso_smiles


RDLogger.logger().setLevel(RDLogger.ERROR)

def get_num_chiral_centers(mol):
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    undef_cc = 0
    for cc in chiral_centers:
        if cc[1] == '?':  
            undef_cc += 1
    return len(chiral_centers), undef_cc
    

def get_num_sp3_centres(mol):
    return sum((x.GetHybridization() == Chem.HybridizationType.SP3) for x in mol.GetAtoms())


def write_error(line, file):
    if file:
        file.write(line + '\n')
        file.flush()
                

def shard(inputs, source, version, output_dir, delimiter,
    name_column=None, skip_lines=0, interval=0, errors_file=None):

    if not os.path.isdir(output_dir):
        utils.log('Creating outdir', output_dir)
        os.mkdir(output_dir)
        
    hacdir = os.path.join(output_dir, source + '_' + version)
    if not os.path.isdir(hacdir):
        utils.log('Creating hacdir', hacdir)
        os.mkdir(hacdir)
    
    count = 0
    errors = 0
    duplicates = 0
    freqs = {}
    files = {}
    
    
    if errors_file:
        err_f = open(os.path.join(output_dir, errors_file), 'w')
    else:
        err_f = None
    
    for input in inputs:
        utils.log_dm_event("Processing", input)

        with (gzip.open(input, 'rt') if input.endswith('.gz') else open(input, 'rt')) as f:
        
            # skip header lines
            if skip_lines:
                for i in range(skip_lines):
                    line = f.readline()
                    
            for line in f:
                if line:
                    count += 1
                    
                    if interval and count % interval == 0:
                        utils.log_dm_event("Processed {} records".format(count))
                    
                    # read the line
                    line = line.strip()
                    tokens = line.split(delimiter)
                    o_smi = tokens[0]
                    
                    # standardize the molecule
                    try:
                        std_smi, mol = standardize_to_iso_smiles(o_smi)
                    except Exception as ex:
                        errors += 1
                        utils.log_dm_event('Error during standardisation of', line)
                        write_error(line, err_f)
                        continue
                    
                    if not mol:
                        errors += 1
                        utils.log_dm_event("Failed to process", line)
                        write_error(line, err_f)
                        continue
                        
                    if name_column is not None:
                        id = tokens[name_column]
                    else:
                        id = str(count)
                        
                    
                    # generate the sha256 digest of the standardized smiles
                    sha256_digest = sha256(bytes(std_smi, 'utf-8')).hexdigest()
                    hac = mol.GetNumHeavyAtoms()
                    if hac > 99:
                        utils.log_dm_event('Molecule has {} atoms - skipping'.format(hac))
                        continue
                            
                    # find if the molecule exists, or if not create the sharded data
                    parts = [output_dir, 'sha256']
                    parts.extend(utils.get_path_from_digest(sha256_digest))
                    path = os.path.join(*parts)
                    if not os.path.isdir(path):
                        os.makedirs(path)
                    molj = os.path.join(path, sha256_digest + '.json')
                    if os.path.exists(molj):
                        # molecule already exists
                        #utils.log("Molecule already exists", id, std_smi)
                        duplicates += 1
                        with open(molj) as f:
                            data = json.load(f)
                  
                        if 'smiles' in data:
                            if data['smiles'] != std_smi:
                                utils.log_dm_event('WARNING: SMILES are inconsistent {} {}'.format(data['smiles'], std_smi))
                                continue
                        else:
                            data['smiles'] = std_smi
                    else:
                        # new molecule - need to create everything
                        data = {}
                        data['smiles'] = std_smi
                        data['uuid'] = str(uuid.uuid4())
                        
                    if not 'suppliers' in data:
                        data['suppliers'] = {}
                    if name_column is not None:
                        if source not in data['suppliers']:
                            data['suppliers'][source] = []
                        if id not in data['suppliers'][source]:
                            data['suppliers'][source].append(id)
                    with open(molj,'w') as f:
                        json.dump(data, f)    
                        
                    uid = data['uuid']   
                    mol.SetProp('_Name', uid)
                    
                    # calculate the molecular props
                    num_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
                    num_rings = rdMolDescriptors.CalcNumRings(mol)
                    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
                    num_cc, num_undef_cc = get_num_chiral_centers(mol)
                    num_sp3 = get_num_sp3_centres(mol)
                    
                    # write the vendor specifc data files (hac shards)
                    if hac in freqs:
                        freqs[hac] = freqs[hac] + 1
                        outfile = files[hac]
                    else:
                        freqs[hac] = 1
                        path = os.path.join(hacdir, str(hac).zfill(2) + '.smi')
                        outfile = open(path, 'w')
                        outfile.write('smiles\tuuid\tid\torig_smiles\tsha256\thac\trot_bonds\tring_count\t' + 
                            'aromatic_ring_count\tchiral_centres\tundefined_chiral_centres\tnum_sp3\n')
                        files[hac] = outfile
                    outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        std_smi, uid, id, o_smi, sha256_digest, hac, num_rot_bonds, num_rings, num_aromatic_rings, num_cc, num_undef_cc, num_sp3))
                        
                    
    
    if err_f:
        err_f.close()
                    
    for f in files.values():
        f.close()
            
    for k in sorted(freqs):
        utils.log(k, freqs[k])
            
    return count, duplicates, errors


def main():

    # Example usage:
    #   ./shard.py -i data/100000.smi -s chemspace -v feb2021 -o testdir -n 1 --interval 10000

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Shard HAC')
    parser.add_argument('-i', '--input', nargs='+', help="Input file(s) as SMILES")
    parser.add_argument('-s', '--source', required=True, help="Molecule source e.g. chemspace")
    parser.add_argument('-v', '--version', required=True, help="Source version e.g. jan_2021")
    parser.add_argument('-o', '--outdir', required=True, help="Dir for the molecules output")
    parser.add_argument('-d', '--delimiter', default='\t', help="Delimiter")
    parser.add_argument('-n', '--name-column', type=int, help="Column for name field (zero based)")
    parser.add_argument('--skip-lines', default=0, type=int, help="Skip this many lines e.g. use 1 for skipping a header line")
    parser.add_argument('-e', '--errors-file', help="Optional file to write bad lines to")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log_dm_event("shard: ", args)
    
    t0 = time.time()
    count, duplicates, errors = shard(args.input, args.source, args.version, args.outdir, args.delimiter, name_column=args.name_column,
        interval=args.interval, skip_lines=args.skip_lines, errors_file=args.errors_file)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    utils.log_dm_event('Processed {} records in {} seconds. {} duplicates, {} errors.'.format(count, duration_s, duplicates, errors))
    

if __name__ == "__main__":
    main()
    
