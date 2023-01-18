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
Synthetic accessibility score.
This is based on the work of Peter Ertl and Greg Landrum that can be found here:
https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score

That in turn is based on this paper:
Peter Ertl and Ansgar Schuffenhauer
Journal of Cheminformatics 1:8 (2009)
http://www.jcheminf.com/content/1/1/8
"""

import argparse, time, os, pickle, traceback, math

import utils, rdkit_utils
from dm_job_utilities.dm_log import DmLog

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

_fscores = None


def readFragmentScores(name='fpscores'):
    import gzip
    global _fscores
    # generate the full path filename:
    if name == "fpscores":
        name = os.path.join(os.path.dirname(__file__), name)
    data = pickle.load(gzip.open('%s.pkl.gz' % name))
    outDict = {}
    for i in data:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _fscores = outDict


def numBridgeheadsAndSpiro(mol, ri=None):
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return nBridgehead, nSpiro


def calculateScore(m):
    if _fscores is None:
        readFragmentScores()

    # fragment score
    fp = rdMolDescriptors.GetMorganFingerprint(m,
                                               2)  # <- 2 is the *radius* of the circular fingerprint
    fps = fp.GetNonzeroElements()
    score1 = 0.
    nf = 0
    for bitId, v in fps.items():
        nf += v
        sfp = bitId
        score1 += _fscores.get(sfp, -4) * v
    score1 /= nf

    # features score
    nAtoms = m.GetNumAtoms()
    nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
    ri = m.GetRingInfo()
    nBridgeheads, nSpiro = numBridgeheadsAndSpiro(m, ri)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgeheads + 1)
    macrocyclePenalty = 0.
    # ---------------------------------------
    # This differs from the paper, which defines:
    #  macrocyclePenalty = math.log10(nMacrocycles+1)
    # This form generates better results when 2 or more macrocycles are present
    if nMacrocycles > 0:
        macrocyclePenalty = math.log10(2)

    score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    # correction for the fingerprint density
    # not in the original publication, added in version 1.1
    # to make highly symmetrical molecules easier to synthetise
    score3 = 0.
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * .5

    sascore = score1 + score2 + score3

    # need to transform "raw" value into scale between 1 and 10
    min = -4.0
    max = 2.5
    sascore = 11. - (sascore - min + 1) / (max - min) * 9.
    # smooth the 10-end
    if sascore > 8.:
        sascore = 8. + math.log(sascore + 1. - 9.)
    if sascore > 10.:
        sascore = 10.0
    elif sascore < 1.:
        sascore = 1.0

    return sascore


def process(input, outfile, delimiter, id_column=None, read_header=False, write_header=False,
            sdf_read_records=100, interval=0):

    utils.expand_path(outfile)

    count = 0
    errors = 0

    # setup the reader
    reader = rdkit_utils.create_reader(input, id_column=id_column, sdf_read_records=sdf_read_records,
                                       read_header=read_header, delimiter=delimiter)
    extra_field_names = reader.get_extra_field_names()

    calc_field_names = ['sa_score']

    # setup the writer
    writer = rdkit_utils.create_writer(outfile, extra_field_names=extra_field_names, calc_prop_names=calc_field_names,
                                       delimiter=delimiter)

    # read the input records and write the output
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break

        mol, smi, id, props = t
        if count == 0 and write_header:
            headers = rdkit_utils.generate_header_values(extra_field_names, len(props), calc_field_names)
            writer.write_header(headers)

        count += 1

        if interval and count % interval == 0:
            DmLog.emit_event("Processed {} records".format(count))
        if count % 50000 == 0:
            # Emit a 'total' cost, replacing all prior costs
            DmLog.emit_cost(count)

        if not mol:
            errors += 1
            DmLog.emit_event("Failed to process record", count)
            continue

        # calculate the molecular props
        try:
            sa_score = calculateScore(mol)

        except:
            errors += 1
            DmLog.emit_event('Failed to process record', count)
            traceback.print_exc()
            continue

        # write the output
        writer.write(smi, mol, id, props, [sa_score])

    writer.close()
    reader.close()

    return count, errors


def main():

    # Example usage:
    #   ./sa_score.py -i data/100.smi -o out.smi

    parser = argparse.ArgumentParser(description='SA Score')
    parser.add_argument('-i', '--input', required=True, help="Input file as SMILES or SDF")
    parser.add_argument('-o', '--outfile', required=True, help="Output file as SMILES or SDF")
    # to pass tab as the delimiter specify it as $'\t' or use one of the symbolic names 'comma', 'tab', 'space' or 'pipe'
    parser.add_argument('-d', '--delimiter', help="Delimiter when using SMILES")
    parser.add_argument('--id-column', help="Column for name field (zero based integer for .smi, text for SDF)")
    parser.add_argument('--read-header', action='store_true',
                        help="Read a header line with the field names when reading .smi or .txt")
    parser.add_argument('--write-header', action='store_true', help='Write a header line when writing .smi or .txt')
    parser.add_argument('--sdf-read-records', default=100, type=int,
                        help="Read this many SDF records to determine field names")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("rdk_props.py: ", args)

    delimiter = utils.read_delimiter(args.delimiter)

    t0 = time.time()
    count, errors = process(args.input, args.outfile, delimiter, id_column=args.id_column,
                            read_header=args.read_header, write_header=args.write_header,
                            sdf_read_records=args.sdf_read_records, interval=args.interval, )
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processed {} records in {} seconds. {} errors.'.format(count, duration_s, errors))
    # Emit final 'total' cost, replacing all prior costs
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
