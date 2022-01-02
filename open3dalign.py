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

"""
3D alignment using Open3DAlign
"""

import argparse, os
import utils
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import rdkit_utils


def process_mol(ref_mol, probe_mol, use_crippen):
    """
    Generate the Open3DAlign alignment and score of a single probe conformer
    :param ref_mol:
    :param probe_mol:
    :param use_crippen:
    :return:
    """
    if use_crippen:
        pyO3A = rdMolAlign.GetCrippenO3A(probe_mol, ref_mol)
    else:
        pyO3A = rdMolAlign.GetO3A(probe_mol, ref_mol)
    if not pyO3A:
        raise ValueError('Failed to read molecule')
    align = pyO3A.Align()
    score = pyO3A.Score()
    return align, score


def process_confs(ref_mol, probe_mol, use_crippen):
    """
    Generate the Open3DAlign alignments and scores for a set of probe conformers
    :param ref_mol:
    :param probe_mol: Can contain multiple conformers
    :param use_crippen:
    :return: tuple of the alignments, tuple of the scores, the best score, the best index
    """
    if use_crippen:
        pyO3As = rdMolAlign.GetCrippenO3AForProbeConfs(probe_mol, ref_mol)
    else:
        pyO3As = rdMolAlign.GetO3AForProbeConfs(probe_mol, ref_mol)

    if len(pyO3As) == 0:
        return None, None, None, None

    best_score = 0
    best_idx = None

    aligns = []
    scores = []
    for idx, pyO3A in enumerate(pyO3As):
        align = pyO3A.Align()
        score = pyO3A.Score()
        aligns.append(align)
        scores.append(score)
        if score > best_score:
            best_score = score
            best_idx = j

    return aligns, scores, best_score, best_idx


def execute(inputs_sdf, query_file, outfile_sdf, use_crippen=False, threshold=None,
            remove_hydrogens=False, interval=None):

    input_count = 0
    output_count = 0
    error_count = 0
    sum_score = 0.0

    q_mol = rdkit_utils.rdk_read_single_mol(query_file)
    q_mol = Chem.AddHs(q_mol, addCoords=True)

    align_perf, score_perf = process_mol(q_mol, q_mol, use_crippen)
    utils.log_dm_event('Query self-alignment score and align:', score_perf , align_perf)

    utils.log_dm_event('Opening', outfile_sdf, 'as output')
    writer = Chem.SDWriter(outfile_sdf)
    try:
        # read the conformers
        suppl = rdkit_utils.rdk_mol_supplier(inputs_sdf)

        # TODO implement the group_by_field functionality
        # mols_in_group = []
        # mols_to_write = None
        # current_group_field_value = None

        # iterate through the conformers and calculate the alignment
        for conf in suppl:
            conf = Chem.AddHs(conf, addCoords=True)
            align, score = process_mol(q_mol, conf, use_crippen)
            sum_score += score
            rel_score = score / score_perf
            input_count += 1

            if not threshold or threshold < rel_score:
                conf.SetDoubleProp('o3da_score', score)
                conf.SetDoubleProp('o3da_score_rel', rel_score)
                conf.SetDoubleProp('o3da_align', align)

                if remove_hydrogens:
                    conf = Chem.RemoveHs(conf)

                writer.write(conf)
                output_count += 1

            if interval and input_count % interval == 0:
                utils.log_dm_event("Processed {} molecules, {} outputs".format(input_count, output_count))

    finally:
        writer.close()

    mean_score = sum_score / input_count
    return input_count, output_count, error_count, mean_score, score_perf
                    

def write_best_mol(writer, mols):
    mols.sort(key=lambda t: t[0])
    writer.write(mols[0][1])


def main():

    # Example:
    #   python3 open3dalign.py -i database.sdf -q query.mol -o results.sdf

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Open3DAlign')
    parser.add_argument('-i', '--inputs', required=True, help="File with molecules to search (SDF)")
    parser.add_argument('-q', '--query', required=True, help="File with the 3D query molecules (SDF or MOL)")
    parser.add_argument('-o', '--outfile', default='o3da-similarity.sdf', help="Output SD file for results")
    parser.add_argument('-c', '--crippen', action='store_true', help='Include Crippen (lipophilicity) considerations')
    parser.add_argument('-t', '--threshold', type=float, help="Score threshold (0 - 1)")
    parser.add_argument('-r', '--remove-hydrogens', action='store_true', help='Remove hydrogens from the outputs')
    # parser.add_argument('-g', '--group-by-field', help="Field name to group records by and report only the best")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log_dm_event("open3dalign.py: ", args)

    input_count, output_count, error_count, mean_score, score_perf = \
        execute(args.inputs, args.query, args.outfile, args.group_by_field,
                use_crippen=args.crippen, threshold=args.threshold, remove_hydrogens=args.remove_hydrogens,
                interval=args.interval)

    tmpl1 = 'Processed {} conformers. Generated {} outputs. {} errors.'
    tmpl2 = 'Perfect score: {} Average score: {}'
    utils.log_dm_event(tmpl1.format(input_count, output_count, error_count))
    utils.log_dm_event(tmpl2.format(score_perf, mean_score))
    
    
if __name__ == "__main__":
    main()
        