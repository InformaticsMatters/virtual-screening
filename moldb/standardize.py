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


import argparse, os, json, gzip, time, logging
import utils
from dm_job_utilities.dm_log import DmLog

from rdkit import Chem, RDLogger
from standardize_molecule import standardize_to_iso_smiles


RDLogger.logger().setLevel(RDLogger.ERROR)


def write_error(line, file):
    if file:
        file.write(line + '\n')
        file.flush()


def standardize(inputs, output_file, delimiter, name_column=None, skip_lines=0, interval=0, write_header=False, errors_file=None):

    count = 0
    errors = 0

    if errors_file:
        utils.expand_path(errors_file)
        err_f = open(errors_file, 'w')
    else:
        err_f = None

    utils.expand_path(output_file)
    with open(output_file, 'w') as outfile:
        if write_header:
            outfile.write('smiles\torig_smiles\tid\n')

        for input in inputs:
            DmLog.emit_event("Processing", input)

            with (gzip.open(input, 'rt') if input.endswith('.gz') else open(input, 'rt')) as f:

                # skip header lines
                if skip_lines:
                    for i in range(skip_lines):
                        line = f.readline()

                for line in f:
                    if line:
                        count += 1

                        if interval and count % interval == 0:
                            DmLog.emit_event("Processed {} records".format(count))
                            if interval * 10 <= count:
                                interval = interval * 2
                        if count % 50000 == 0:
                            DmLog.emit_cost(count)

                        # read the line
                        line = line.strip()
                        tokens = line.split(delimiter)
                        o_smi = tokens[0]
                        id = tokens[name_column]

                        # standardize the molecule
                        try:
                            std_smi, mol = standardize_to_iso_smiles(o_smi)
                        except Exception as ex:
                            errors += 1
                            DmLog.emit_event('Error during standardisation of', line)
                            write_error(line, err_f)
                            continue

                        if not mol:
                            errors += 1
                            DmLog.emit_event("Failed to process", line)
                            write_error(line, err_f)
                            continue

                        outfile.write("{}\t{}\t{}\n".format(std_smi, o_smi, id))

    if err_f:
        err_f.close()

    return count, errors


def main():

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Standardize')
    parser.add_argument('-i', '--input', nargs='+', help="Input file(s) as SMILES")
    parser.add_argument('-o', '--outfile', default='output.smi', help="Output file")
    parser.add_argument('-d', '--delimiter', default='\t', help="Delimiter")
    parser.add_argument('-n', '--name-column', required=True, type=int, help="Column for name field (zero based)")
    parser.add_argument('--skip-lines', default=0, type=int, help="Skip this many lines e.g. use 1 to skip header line")
    parser.add_argument('--write-header', action='store_true', help="Write a header line")
    parser.add_argument('-e', '--errors-file', help="Optional file to write bad lines to")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("standardize: ", args)

    delimiter = utils.read_delimiter(args.delimiter)

    t0 = time.time()
    count, errors = standardize(args.input, args.outfile, delimiter, name_column=args.name_column,
                                write_header=args.write_header, interval=args.interval, skip_lines=args.skip_lines,
                                errors_file=args.errors_file)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processed {} records in {} seconds., {} errors.'.format(count, duration_s, errors))
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()

