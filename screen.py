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

import argparse, time
import traceback

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols

import utils, rdkit_utils

"""
Filter a file of SMILES using fingerprint based similarity.
One or more query molecules can be specified as SMILES, or the queries can be contained in a SD file or SMILES file. 
Similarities for each query are calculated.
When more than one query is specified the min similarity, max similarity, the arithmetic mean of the similarities, the
geometric mean of the similarities and the product of the similarities are also generated and added to the list of 
similarities.
The inputs are filtered using a threshold score and the index of the similarity score to use.
For instance when specifying the query as SMILES and where there is only a single query molecule there is only a single 
score, so the index to use is 0.
When there are two query molecules specified as SMILES there are two scores plus the extra four. So if you want to 
filter by the sum of the individual scores then use an index of 4. For two inputs the indices are:
- 0 similarity to first molecule
- 1 similarity to second molecule
- 2 minimum of the individual scores
- 3 maximum of the individual scores
- 4 arithmetic mean of the individual scores
- 5 geometric mean of the individual scores
- 6 product of the individual scores

When the queries are specified in a file it is assumed that there will be a relatively large number of query molecules
so the individual scores are not output, nor is the product as that is almost always close to zero. In this case the 
index to use to filter will be:
- 0 minimum of the individual scores
- 1 maximum of the individual scores
- 2 arithmetic mean of the individual scores
- 3 geometric mean of the individual scores

The descriptor and metric to use can be specified. See the descriptors and metrics properties for the permitted values.
When using Morgan fingerprints you need to be careful about which metric you use. By default Morgan fingerprints use
counts, but these can only be used with tanimoto, dice or tversky metrics. If using another metric you must generate
Morgan bit vectors which you do by specifying the nbits parameter e.g. with a value of 1024. 
See http://rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints

When using the tversky metric you can also specify the alpha and beta parameters. Default values are 1 and 0. These
parameters are ignored when using other metrics.  

The output is that same as the input with extra similarity scores appended to each line, with only lines passing the 
threshold filter being written.
"""

descriptors = {
    'maccs':   lambda m, nBits: MACCSkeys.GenMACCSKeys(m),
    'morgan2': lambda m, nBits: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=nBits) if nBits else AllChem.GetMorganFingerprint(m, 2),
    'morgan3': lambda m, nBits: AllChem.GetMorganFingerprintAsBitVect(m, 3, nBits=nBits) if nBits else AllChem.GetMorganFingerprint(m, 3),
    'rdkit':   lambda m, nBits: FingerprintMols.FingerprintMol(m),
}

metrics = {
    'asymmetric': lambda m1, m2, a, b: DataStructs.AsymmetricSimilarity(m1, m2),
    'braunblanquet': lambda m1, m2, a, b: DataStructs.BraunBlanquetSimilarity(m1, m2),
    'cosine': lambda m1, m2, a, b: DataStructs.CosineSimilarity(m1, m2),
    'dice': lambda m1, m2, a, b: DataStructs.DiceSimilarity(m1, m2),
    'kulczynski': lambda m1, m2, a, b: DataStructs.KulczynskiSimilarity(m1, m2),
    'mcconnaughey': lambda m1, m2, a, b: DataStructs.McConnaugheySimilarity(m1, m2),
    'rogotgoldberg': lambda m1, m2, a, b: DataStructs.RogotGoldbergSimilarity(m1, m2),
    'russel': lambda m1, m2, a, b: DataStructs.RusselSimilarity(m1, m2),
    'sokal': lambda m1, m2, a, b: DataStructs.SokalSimilarity(m1, m2),
    'tanimoto': lambda m1, m2, a, b: DataStructs.TanimotoSimilarity(m1, m2),
    'tversky': lambda m1, m2, a, b: DataStructs.TverskySimilarity(m2, m1, a, b)
}


def multiply_list(scores) :
    result = 1
    for x in scores:
        result = result * x
    return result


def calc_geometric_mean(scores):
    total = 1.0
    for score in scores:
        total = total * score
    result = total ** (1.0/len(scores))
    return result


def validate_params(descriptor, metric, alpha, beta, nbits):
    if descriptor == 'morgan2' or descriptor == 'morgan3':
        if not (metric == 'tamimoto' or metric == 'dice' or metric == 'tversky'):
            if not nbits:
                utils.log_dm_event('When using', descriptor, 'descriptor and', metric, 'metric',
                                   'the nbits parameter must be defined')
                exit(1)
            else:
                if metric == 'tversky':
                    utils.log_dm_event('Using {} bit vector and {} metric with alpha={}, beta={}'
                                       .format(descriptor, metric, alpha, beta))
                else:
                    utils.log_dm_event('Using {} bit vector and {} metric'.format(descriptor, metric))
        else:
            if metric == 'tversky':
                utils.log_dm_event('Using {} counts and {} metric with alpha={}, beta={}'
                                   .format(descriptor, metric, alpha, beta))
            else:
                utils.log_dm_event('Using {} counts and {} metric'.format(descriptor, metric))
    elif metric == 'tversky':
        utils.log_dm_event('Using {} descriptor and {} metric with alpha={}, beta={}'
                           .format(descriptor, metric, alpha, beta))
    else:
        utils.log_dm_event('Using {} descriptor and {} metric'.format(descriptor, metric))


def execute(query_smis, query_file, inputfile, outputfile, descriptor, metric,
            delimiter='\t', threshold=0.7, sim_idx=0, read_header=False, write_header=False,
            queries_read_header=False, queries_delimiter=None, alpha=1.0, beta=0.0, nbits=None, interval=None):

    if query_smis and query_file:
        raise ValueError("Specify queries as SMILES or a file, not both")

    validate_params(descriptor, metric, alpha, beta, nbits)

    count = 0
    hits = 0
    errors = 0

    my_descriptor = descriptors[descriptor]
    my_metric = metrics[metric]

    q_mols = []
    if query_smis:
        for smi in query_smis:
            mol = Chem.MolFromSmiles(smi)
            if not mol:
                raise ValueError('Failed to read query smiles:', smi)
            q_mols.append(mol)
    elif query_file:
        reader = rdkit_utils.create_reader(query_file, read_header=queries_read_header, delimiter=queries_delimiter)
        q_count = 0
        while True:
            q_count += 1
            t = reader.read()
            # break if no more data to read
            if not t:
                break
            mol, smi, id, props = t
            if not mol:
                raise ValueError('Failed to read query molecule:', q_count)
            q_mols.append(mol)
        utils.log_dm_event('Read', len(q_mols), 'query molecules')

    q_fps = []
    for q_mol in q_mols:
        q_fp = my_descriptor(q_mol, nbits)
        q_fps.append(q_fp)

    with open(outputfile, 'wt') as outf:
        with open(inputfile) as inf:
            if read_header:
                line = next(inf)
                line = line.strip()
                headers = line.split(delimiter)
            else:
                headers = None

            for line in inf:
                count += 1

                if interval and count % interval == 0:
                    utils.log_dm_event("Processed {} records, {} hits".format(count, hits))

                line = line.strip()
                tokens = line.split(delimiter)
                smi = tokens[0]
                if write_header and count == 1:
                    if headers is None:
                        # we need to generate generic field names
                        headers = []
                        for i, t in enumerate(tokens):
                            if i == 0:
                                headers.append('smiles')
                            else:
                                headers.append('field' + str(i + 1))
                    if query_smis:
                        for i in range(len(query_smis)):
                            headers.append('score_' + str(i + 1))
                    if query_smis and len(query_smis) > 1:
                        headers.extend(['score_min', 'score_max', 'score_amean', 'score_gmean', 'score_prod'])
                    if query_file:
                        headers.extend(['score_min', 'score_max', 'score_amean', 'score_gmean'])

                    outf.write(delimiter.join(headers) + '\n')
                try:
                    t_mol = Chem.MolFromSmiles(smi)
                    if not t_mol:
                        errors += 1
                        utils.log_dm_event('Failed to process molecule', count, smi)
                        continue
                    t_fp = my_descriptor(t_mol, nbits)
                    sims = []
                    for q_fp in q_fps:
                        sims.append(my_metric(q_fp, t_fp, alpha, beta))

                    if len(sims) > 1:
                        min_sims = min(sims)
                        max_sims = max(sims)
                        amean_sims = sum(sims) / len(q_mols)
                        gmean_sims = calc_geometric_mean(sims)
                        prod_sims = multiply_list(sims)
                        if query_file:
                            # clear the sims as we don't want to write individual scores when using a query file
                            sims = []
                        sims.append(min_sims)
                        sims.append(max_sims)
                        sims.append(amean_sims)
                        sims.append(gmean_sims)
                        if query_smis:
                            sims.append(prod_sims)

                    sim = sims[sim_idx]
                    if sim > threshold:
                        hits += 1
                        #utils.log(smi, sims)
                        outl = line
                        for sim in sims:
                            outl += delimiter + str(sim)

                        outf.write(outl + '\n')
                except Exception as e:
                    errors += 1
                    utils.log('Failed to process molecule', count, smi)
                    traceback.print_exc()

    return count, hits, errors


### start main execution #########################################

def main():

    # Example:
    #   python3 screen.py --smiles 'O=C(Nc1ccc(Cl)cc1)c1ccccn1' --input data/10000.smi --delimiter tab -o foo.smi\
    #     -d morgan2 -m tanimoto

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='screen')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('-s', '--smiles', nargs='+', help="Query SMILES")
    inputs.add_argument('-q', '--queries', help="File with query molecules")
    parser.add_argument('--queries-delimiter', help="Delimiter for queries file (text format)")
    parser.add_argument('--queries-read-header', action='store_true',
                        help="Does the queries file contain a header line (text format)")
    parser.add_argument('-i', '--input', required=True, help="SMILES file with molecules to search")
    parser.add_argument('--delimiter', default='\t', help="Delimiter")
    parser.add_argument('--read-header', action='store_true', help="Read a header line with the field names")
    parser.add_argument('-o', '--output', required=True, help="Output file as SMILES")
    parser.add_argument('--write-header', action='store_true', help='Write a header line')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(descriptors.keys()), default='rdkit',
                        help='Descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric', type=str.lower, choices=list(metrics.keys()), default='tanimoto',
                        help='Similarity metric (default tanimoto)')
    parser.add_argument("--threshold", type=float, default=0.7, help="Similarity threshold")
    parser.add_argument("--sim-index", type=int, default=0, help="Similarity score index")
    parser.add_argument("--alpha", type=float, default=1.0, help="Tversky alpha parameter")
    parser.add_argument("--beta", type=float, default=0.0, help="Tversky beta parameter")
    parser.add_argument("--nbits", type=int, default=None, help="Number of bits if using Morgan  as bit vector e.g. 1024")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log_dm_event("screen: ", args)

    delimiter = utils.read_delimiter(args.delimiter)
    queries_delimiter = utils.read_delimiter(args.queries_delimiter)

    start = time.time()
    input_count, hit_count, error_count = \
        execute(args.smiles, args.queries, args.input, args.output, args.descriptor, args.metric,
                threshold=args.threshold, sim_idx=args.sim_index, delimiter=delimiter,
                read_header=args.read_header, write_header=args.write_header,
                queries_read_header=args.queries_read_header, queries_delimiter=queries_delimiter,
                alpha=args.alpha, beta=args.beta, nbits=args.nbits, interval=args.interval
                )
    end = time.time()
    duration_s = int(end - start)
    if duration_s < 1:
        duration_s = 1

    utils.log_dm_event(input_count, 'inputs,', hit_count, 'hits,', error_count, 'errors.',
                         'Time (s):', duration_s)


if __name__ == "__main__":
    main()
