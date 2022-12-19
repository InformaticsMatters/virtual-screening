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

import utils
from dm_job_utilities.dm_log import DmLog
from openbabel import pybel, OBConformerSearch
from jinja2 import Template
import argparse, os, time, tempfile, glob
import subprocess

template = Template("""
score donor donor -6.00 0.25
score acceptor acceptor -6.00 0.25
score donacc donacc -6.00 0.25
score donacc acceptor -3.00 0.25
score donacc donor -3.00 0.25
score nonpolar_noring nonpolar_noring -0.25 0.25

ringscore -10.00 0.25

alignment_torsion_weight {{ alignment_torsion_weight }}

aco_sigma 1.0
aco_ants 20
aco_evap 0.2

pher_smoothings 3
pher_pbest 0.5
pher_force_descent_update 5
simplex_tolerance 0.02
refine_simplex_tolerance 0.0001

{% for frag in fragments %}
ligand_file {{ frag }} fixed{% endfor %}
ligand_file {{ molecule }}

cluster_structures {{ cluster_structures }}
cluster_rmsd {{ cluster_rmsd }}

output_dir {{ output_dir }}

""")


def process(inputs, fragments, outputfile, alignment_torsion_weight=20, cluster_structures=5, cluster_rmsd=2,
            threshold=None, gen_coords=False, dir=None, interval=None):

    if dir:
        if not os.path.exists(dir):
            os.makedirs(dir, exist_ok=True)
        num_mols, num_written, num_errors = _do_processing(dir, inputs, fragments, outputfile, gen_coords,
                                              alignment_torsion_weight, cluster_structures,
                                              cluster_rmsd, threshold, interval)
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            utils.log('created temporary directory', tmpdirname)
            num_mols, num_written, num_errors = _do_processing(tmpdirname, inputs, fragments, outputfile, gen_coords,
                                                  alignment_torsion_weight, cluster_structures,
                                                  cluster_rmsd, threshold, interval)

    return num_mols, num_written, num_errors

def _write_frag_mol2(mol, dir, index):
    mol.addh()
    p = os.path.join(dir, 'frag' + str(index) + '.mol2')
    mol.write(format='mol2', filename=p)
    return p


def _do_processing(dir, inputs, fragments, outputfile, gen_coords, alignment_torsion_weight, cluster_structures,
                   cluster_rmsd, threshold, interval):
    num_mols = 0
    num_frags = 0
    num_written = 0
    num_errors = 0

    frag_paths = []
    mol_paths = []

    # write the fragments as mol2 files
    for frag in fragments:
        if frag.endswith('.mol'):
            mol = next(pybel.readfile("mol", frag))
            p = _write_frag_mol2(mol, dir, num_frags)
            frag_paths.append(p)
            num_frags += 1
        elif frag.endswith('.sdf'):
            mols = pybel.readfile("sdf", frag)
            for mol in mols:
                p = _write_frag_mol2(mol, dir, num_frags)
                frag_paths.append(p)
                num_frags += 1
        else:
            raise ValueError('Fragments must be .mol or .sdf. Found', frag)

    utils.log('Fragments:', frag_paths)
    DmLog.emit_event("Found {} fragments".format(num_frags))

    # write the inputs as mol2 files
    for mol in inputs:
        if gen_coords:
            mol.make3D()
        props = {}
        props.update(mol.data)
        if 'MOL Chiral Flag' in props:
            del props['MOL Chiral Flag']
        mol.addh()
        p = os.path.join(dir, 'mol' + str(num_mols) + '.mol2')
        mol.write(format='mol2', filename=p)
        chg_block = _read_charge_block(p)
        mol_paths.append([p, props, chg_block])
        num_mols += 1

        if interval and num_mols % interval == 0:
            DmLog.emit_event("Processed {} records".format(num_mols))

    # run pharmACOphore
    for i, mol_path in enumerate(mol_paths):
        # create the config file
        path = mol_path[0]
        output_dir = os.path.join(dir, 'results' + str(i))
        mol_path.append(output_dir)
        content = template.render(alignment_torsion_weight=alignment_torsion_weight,
                                  cluster_structures=cluster_structures,
                                  cluster_rmsd=cluster_rmsd,
                                  fragments=frag_paths,
                                  molecule=path,
                                  output_dir=output_dir)

        config = os.path.join(dir, 'align' + str(i) + '.config')
        with open(config, 'wt') as writer:
            writer.write(content)

        # run plants
        cmd = ['plants', '--mode', 'align', config]
        #utils.log("CMD: " + " ".join(cmd))
        proc = subprocess.run(cmd, capture_output=True)

        utils.log('Aligned molecule', i)

    # collate the results to a SD file
    utils.expand_path(outputfile)
    utils.log(('Writing output to', outputfile))
    with pybel.Outputfile('sdf', outputfile, overwrite=True) as writer:
        for mol_path in mol_paths:
            props = mol_path[1]
            chg_block = mol_path[2]
            output_dir = mol_path[3]
            utils.log('Handling', output_dir)

            # read the scores
            csv_file = os.path.join(output_dir, 'ranking.csv')
            if not os.path.exists(csv_file):
                num_errors += 1
            else:
                scores = []
                with open(csv_file) as reader:
                    header = reader.readline()
                    while True:
                        line = reader.readline()
                        if not line:
                            break
                        else:
                            tokens = line.split(',')
                            scores.append(tokens[1])

                # read the mol2 files and write to SDF
                for i, f in enumerate(glob.glob(os.path.join(output_dir, 'aligned_mol*_conf_*.mol2'))):
                    #utils.log('Handling molecule', i, f)
                    if not threshold or float(scores[i]) <= threshold:
                        mol2block = _write_charge_block(f, chg_block)
                        mol = pybel.readstring("mol2", mol2block)
                        mol.removeh()
                        mol.data.update(props)
                        mol.data['PH4_SCORE'] = scores[i]
                        writer.write(mol)
                        num_written += 1

    return num_mols, num_written, num_errors


def _write_charge_block(mol2file, chg_block):
    lines = []
    with open(mol2file, 'rt') as reader:
        while True:
            line = reader.readline()
            if not line:
                break
            elif line.strip() == '@<TRIPOS>BOND':
                lines.extend(chg_block)
                lines.append(line)
            elif line.strip() == 'NO_CHARGES':
                lines.append('GASTEIGER\n')
            else:
                lines.append(line)
    return ''.join(lines)


def _read_charge_block(mol2file):
    started = False
    block = []
    with open(mol2file, 'rt') as reader:
        while True:
            line = reader.readline()
            if not line:
                break
            else:
                if started:
                    if line[0] == '@':
                        break
                    else:
                        block.append(line)
                else:
                    if line.strip() == '@<TRIPOS>UNITY_ATOM_ATTR':
                        block.append(line)
                        started = True

    return block


def smiles_generator(smiles):
    for smi in smiles:
        mol = pybel.readstring('smi', smi)
        yield mol


def file_generator(filename):
    ftype = 'sdf' if filename.endswith('.sdf') else 'smi'
    supplr = pybel.readfile(ftype, filename)
    return supplr


def main():

    # Example usages:
    #   ./pharmacophore.py -i data/candidates.sdf -f data/Mpro-x0107_0A.mol data/Mpro-x1382_0A.mol --outfile out.sdf
    #   ./pharmacophore.py -i data/candidates.sdf -f data/fragments.sdf --outfile out.sdf
    #   python /code/pharmacophore.py -i /data/mpro-merge-enumerated.sdf -f /data/data/Mpro-x0107_0A.mol /data/data/Mpro-x1382_0A.mol --outfile out.sdf --interval 5000 -c 1

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='3D flexible alignment using pharmACOphore')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('-s', '--smiles', nargs='+', help="Input SMILES")
    inputs.add_argument('-i', '--input', help="Input file (.sdf, .smi)")
    parser.add_argument('-f', '--fragments', nargs='+', required=True, help="Molfiles or SD files with fragments")
    parser.add_argument('-o', '--outfile', required=True, help="Output file SDF")
    parser.add_argument('-w', '--work-dir', help="Directory to work in. If not defined a temp dir is used")
    parser.add_argument('-g', '--gen-coords', action='store_true', help="Generate 3D coordinates for inputs")
    parser.add_argument('-t', '--torsion-weight', type=float, default=20,
                        help="Weight for torsional energy contribution")
    parser.add_argument('--threshold', type=float, help="Score threshold")
    parser.add_argument('-r', '--rmsd', type=float, default=2, help="RMSD pruning threshold")
    parser.add_argument('-c', '--count', type=int, default=10, help="Number of alignments to retain")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("pharmacophore.py: ", args)

    if args.smiles:
        supplr = smiles_generator(args.smiles)
    elif args.input:
        supplr = file_generator(args.input)

    t0 = time.time()
    count, written, errors = process(supplr, args.fragments, args.outfile, dir=args.work_dir,
                            cluster_structures=args.count, cluster_rmsd=args.rmsd, gen_coords=args.gen_coords,
                            alignment_torsion_weight=args.torsion_weight, threshold=args.threshold,
                            interval=args.interval)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processed {} records in {} seconds. {} results, {} errors.'.format(count, duration_s, written, errors))
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
