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
Determine shape similarity using the USR family of tools:
- USR: Ballester PJ, Richards WG (2007). DOI: 10.1002/jcc.20681
- Electroshape: Armstrong, M. S. et al. (2010). DOI: 10.1007/s10822-010-9374-0
- USRCAT: Schreyer A, Blundell T (2012). DOI: 10.1186/1758-2946-4-27

Inputs are a query molecule (3D) and a 3D SD file containing conformers to be screened.
Whilst any 3D SD file can in principle be used, this tool is optimized to use the tooling in this project that can be
used to generate the conformers using our standard conventions. See le_conformers.py (and the le_conformers.nf workflow)
and the assemble_conformers.py for more details on how to conform to these conventions.

The output is the input molecules that pass the shape similarity threshold. The similarity score is added as an
additional field. The most similar record among a consecutive set of records, identified by the group-by-field parameter,
can, and usually should, be output. If conformers are generated using the above mentioned tooling the value for the
group-by-field parameter should typically be std_smi.
"""

import argparse, time
import utils
from dm_job_utilities.dm_log import DmLog

from oddt import toolkit
from oddt import shape


def execute(inputs_sdf, queries_file, outfile_sdf, method, group_by_field, threshold, interval=None):

    utils.expand_path(outfile_sdf)

    output_count = 0
    error_count = 0
    sum_similarity = 0.0
    input_count = 0
    hit_count = 0

    if method == 'usr':
        method_func = shape.usr
    elif method == 'electroshape':
        method_func = shape.electroshape
    elif method == 'usrcat':
        method_func = shape.usr_cat
    else:
        raise ValueError('Unsupported method: ' + method)

    # read the query molecule. Can be SDF or Mol format
    # if SDF then the first molecule is used.
    qmol = next(toolkit.readfile('sdf', queries_file))
    qmol.removeh()
    qshape = method_func(qmol)

    DmLog.emit_event('Opening', outfile_sdf, 'as output')
    writer = toolkit.Outputfile('sdf', outfile_sdf, overwrite=True)
    try:
        # read the conformers
        confs = toolkit.readfile('sdf', inputs_sdf)

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
            input_count += 1

            if interval and input_count % interval == 0:
                DmLog.emit_event("Processed {} molecules, {} hits. {} outputs".format(input_count, hit_count, output_count))
            if input_count % 10000 == 0:
                DmLog.emit_cost(input_count)

            if similarity > threshold:
                hit_count += 1
                conf.data[method + '_similarity'] = similarity
                # if 'enum_smi' in conf.data:
                #     utils.log(input_count, similarity, conf.data['enum_smi'])
                # else:
                #     utils.log(input_count, similarity)
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

    mean_similarity = sum_similarity / input_count
    return input_count, output_count, error_count, mean_similarity
                    

def write_best_mol(writer, mols):
    mols.sort(key=lambda t: t[0])
    writer.write(mols[0][1])


def main():

    # Example:
    #   python3 usr.py -i database.sdf -q query.mol -o results.sdf -t 0.6 -m usr

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Ultrafast Shape Recognition')
    parser.add_argument('-i', '--inputs', required=True, help="File with molecules to search")
    parser.add_argument('-q', '--query', required=True, help="File with the 3D query molecules (SDF or MOL)")
    parser.add_argument('-o', '--outfile', default='usr-similarity.sdf', help="Output SD file for results")
    parser.add_argument('-t', '--threshold', required=True, type=float, help="Score threshold")
    parser.add_argument('-m', '--method', required=True, choices=['usr', 'electroshape', 'usrcat'],
                        help='Shape method [usr, electroshape, usrcat]')
    parser.add_argument('-g', '--group-by-field', help="Field name to group records by and report only the best")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("usr.py: ", args)

    start = time.time()
    input_count, output_count, error_count, mean_similarity = \
        execute(args.inputs, args.query, args.outfile, args.method, args.group_by_field,
                args.threshold, interval=args.interval)
    end = time.time()
    duration_s = int(end - start)
    if duration_s < 1:
        duration_s = 1

    tmpl = 'Processed {} conformers. Generated {} outputs. {} errors. Average similarity is {}. Time (s): {}'
    DmLog.emit_event(tmpl.format(input_count, output_count, error_count, mean_similarity, duration_s))
    DmLog.emit_cost(input_count)
    
    
if __name__ == "__main__":
    main()
