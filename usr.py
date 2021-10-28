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


import argparse, os
import utils
from utils import get_path_from_digest

from oddt import toolkit
from oddt import shape


def execute(inputs_smi, queries_file, outfile_sdf, data_dir, method, group_by_field, threshold, interval=None):

    input_count = 0
    output_count = 0
    error_count = 0
    sum_similarity = 0.0
    similarity_count = 0

    if method == 'usr':
        method_func = shape.usr
    elif method == 'electroshape':
        method_func = shape.electroshape
    elif method == 'usrcat':
        method_func = shape.usr_cat
    else:
        raise ValueError('Unsupported method: ' + method)
    
    utils.log_dm_event('Processing file', inputs_smi)

    # Read the inputs and cache as a list of tuples
    inputs = []
    with open(inputs_smi) as inf:
        for line in inf:
            input_count += 1
            tokens = line.strip().split('\t')
            smi = tokens[0]
            uid = tokens[1]
            digest = tokens[2]
            inputs.append((smi, uid, digest))

    utils.log_dm_event('Read', input_count, 'inputs')

    # read the query molecule. Can be SDF or Mol format
    # if SDF then the first molecule is used.
    qmol = next(toolkit.readfile('sdf', queries_file))
    qmol.removeh()
    qshape = method_func(qmol)

    utils.log_dm_event('Opening', outfile_sdf, 'as output')
    writer = toolkit.Outputfile('sdf', outfile_sdf, overwrite=True)
    try:
        # iterate through the inputs
        for input in inputs:
            parts = [data_dir]
            digest = input[2]
            parts.extend(get_path_from_digest(digest))
            path = os.path.join(*parts)
            if not os.path.isdir(path):
                utils.log_dm_event('WARNING, path', path, 'not found')
                error_count += 1
                continue

            # generate the conformers file name and check it exists
            confs_sdf = os.path.join(path, digest + '_le_confs.sdf')
            if not os.path.isfile(confs_sdf):
                utils.log_dm_event('WARNING, path', confs_sdf, 'not found')
                error_count += 1
                continue

            # read the conformers as a list
            confs = list(toolkit.readfile('sdf', confs_sdf))

            # iterate through the conformers and calculate the similarity
            mols_in_group = []
            mols_to_write = None
            current_group_field_value = None

            for conf in confs:
                if mols_to_write:
                    write_best_mol(writer, mols_to_write)
                    mols_to_write = None
                    output_count += 1

                conf.removeh()
                cshape = method_func(conf)
                similarity = shape.usr_similarity(qshape, cshape)
                sum_similarity += similarity
                similarity_count += 1

                if interval and similarity_count % interval == 0:
                    utils.log_dm_event("Processed {} molecules".format(similarity_count))

                if similarity > threshold:
                    conf.data[method + '_similarity'] = similarity
                    utils.log(similarity, input[0])
                    if group_by_field:
                        if group_by_field in conf.data:
                            value = conf.data[group_by_field]
                            if current_group_field_value and current_group_field_value != value:
                                mols_to_write = mols_in_group.copy()
                                mols_in_group = [(similarity, conf)]
                            else:
                                mols_in_group.append((similarity, conf))
                            current_group_field_value = value
                        else:
                            raise ValueError('Grouping field', group_by_field, 'not present')
                    else:
                        mols_to_write = [(similarity, conf)]

            if mols_in_group:
                write_best_mol(writer, mols_in_group)
                output_count += 1

    finally:
        writer.close()

    mean_similarity = sum_similarity / similarity_count
    return input_count, similarity_count, output_count, error_count, mean_similarity
                    

def write_best_mol(writer, mols):
    mols.sort(key=lambda t: t[0])
    writer.write(mols[0][1])


def main():

    # Example:
    #   python3 usr.py -i foo.smi -q foo.mol -t 0.6 -m usr

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Ultrafast Shape Recognition')
    parser.add_argument('-i', '--inputs', required=True, help="File with inputs")
    parser.add_argument('-q', '--query', required=True, help="File with the 3D query molecules (SDF or MOL)")
    parser.add_argument('--outfile', default='usr-similarity.sdf', help="Output SD file for results")
    parser.add_argument('-d', '--data-dir', default='molecules/sha256', help="Directory with data")
    parser.add_argument('-t', '--threshold', required=True, type=float, help="Score threshold")
    parser.add_argument('-m', '--method', required=True, choices=['usr', 'electroshape', 'usrcat'],
                        help='Shape method [usr, electroshape, usrcat]')
    parser.add_argument('-g', '--group-by-field', help="Field name to group records by and report only the best")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log("usr.py: ", args)

    input_count, similarity_count, output_count, error_count, mean_similarity = \
        execute(args.inputs, args.queries, args.outfile, args.data_dir, args.method, args.group_by_field,
                args.threshold, interval=args.interval)

    tmpl = 'Processed {} inputs and {} conformers. Generated {} outputs. {} errors. Average similarity is {}'
    utils.log_dm_event(tmpl.format(
        input_count, similarity_count, output_count, error_count, mean_similarity))
    
    
if __name__ == "__main__":
    main()
        
