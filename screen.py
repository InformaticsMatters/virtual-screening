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

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols

import utils

"""
Filter a file of SMILES using fingerprint based similarity.
One or more query molecules can be specified as SMILES.
Similarities for each query are calculated.
When more than one query is specified the min similarity, max similarity, the sum of the similarities and the product
of the similarities are also generated and added to the list of similarities.
The inputs are filtered using a threshold score and the index of the similarity score to use.
For instance when there is only a single query molecule there is only a single score, so the index to use is 0.
When there are two query molecules there are two scores plus the extra four. So if you want to filter by the sum of the
individual scores then use an index of 4. For two inputs the indices are:
- 0 similarity to first molecule
- 1 similarity to second molecule
- 2 minimum of the individual scores
- 3 maximum of the individual scores
- 4 arithmetic mean of the individual scores
- 5 geometric mean of the individual scores
- 6 product of the individual scores

The descriptor and metric to use can be specified. See the descriptors and metrics properties for the permitted values.

The output is that same as the input with extra similarity scores appended to each line, with only lines passing the 
threshold filter being written.
"""

descriptors = {
    'maccs':   lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan2': lambda m: AllChem.GetMorganFingerprint(m, 2),
    'morgan3': lambda m: AllChem.GetMorganFingerprint(m, 3),
    'rdkit':   lambda m: FingerprintMols.FingerprintMol(m),
}

metrics = {
    'asymmetric': DataStructs.AsymmetricSimilarity,
    'braunblanquet': DataStructs.BraunBlanquetSimilarity,
    'cosine': DataStructs.CosineSimilarity,
    'dice': DataStructs.DiceSimilarity,
    'kulczynski': DataStructs.KulczynskiSimilarity,
    'mcconnaughey': DataStructs.McConnaugheySimilarity,
    'rogotgoldberg': DataStructs.RogotGoldbergSimilarity,
    'russel': DataStructs.RusselSimilarity,
    'sokal': DataStructs.SokalSimilarity,
    'tanimoto': DataStructs.TanimotoSimilarity,
    'tversky': lambda m1, m2: DataStructs.TverskySimilarity(m1, m2, 1, 0),
    'tversky_inverse': lambda m1, m2: DataStructs.TverskySimilarity(m2, m1, 1, 0)
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


def execute(query_smis, inputfile, outputfile, descriptor, metric,
            delimiter='\t', threshold=0.7, sim_idx=0,
            read_header=False, write_header=False, interval=None):

    count = 0
    hits = 0
    errors = 0

    descriptor = descriptors[descriptor]
    metric = metrics[metric]

    q_fps = []
    for q_smi in query_smis:
        q_mol = Chem.MolFromSmiles(q_smi)
        if not q_mol:
            raise ValueError('Failed to read query SMILES')
        q_fp = descriptor(q_mol)
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
                    for i in range(len(query_smis)):
                        headers.append('score_' + str(i + 1))
                    if len(query_smis) > 1:
                        headers.extend(['score_min', 'score_max', 'score_amean', 'score_gmean', 'score_prod'])
                    outf.write(delimiter.join(headers) + '\n')
                try:
                    t_mol = Chem.MolFromSmiles(smi)
                    if not t_mol:
                        errors += 1
                        utils.log_dm_event('Failed to process molecule', count, smi)
                        continue
                    t_fp = descriptor(t_mol)
                    sims = []
                    for q_fp in q_fps:
                        sims.append(metric(q_fp, t_fp))

                    if len(sims) > 1:
                        min_sims = min(sims)
                        max_sims = max(sims)
                        amean_sims = sum(sims) / len(query_smis)
                        gmean_sims = calc_geometric_mean(sims)
                        prod_sims = multiply_list(sims)
                        sims.append(min_sims)
                        sims.append(max_sims)
                        sims.append(amean_sims)
                        sims.append(gmean_sims)
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

    return count, hits, errors



### start main execution #########################################

def main():

    # Example:
    #   python3 screen.py --smiles 'O=C(Nc1ccc(Cl)cc1)c1ccccn1' --input molecules.smi --delimiter tab\
    #     -d morgan2 -m tanimoto

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='screen')
    parser.add_argument('-s', '--smiles', nargs='+', required=True, help="Query SMILES")
    parser.add_argument('-i', '--input', required=True, help="SMILES file with targets")
    parser.add_argument('--delimiter', default='\t', help="Delimiter")
    parser.add_argument('-o', '--output', required=True, help="Output file as SMILES")
    parser.add_argument('--read-header', action='store_true', help="Read a header line with the field names")
    parser.add_argument('--write-header', action='store_true', help='Write a header line')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(descriptors.keys()), default='rdkit', help='Descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric', type=str.lower, choices=list(metrics.keys()), default='tanimoto', help='Similarity metric (default tanimoto)')
    parser.add_argument("--threshold", type=float, default=0.7, help="Similarity threshold")
    parser.add_argument("--sim-index", type=int, default=0, help="Similarity score index")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    utils.log_dm_event("screen: ", args)

    delimiter = utils.read_delimiter(args.delimiter)

    start = time.time()
    input_count, hit_count, error_count = \
        execute(args.smiles, args.input, args.output, args.descriptor, args.metric,
                threshold=args.threshold, sim_idx=args.sim_index, delimiter=delimiter,
                read_header=args.read_header, write_header=args.write_header, interval=args.interval
                )
    end = time.time()
    duration_s = int(end - start)
    if duration_s < 1:
        duration_s = 1

    utils.log_dm_event(input_count, 'inputs,', hit_count, 'hits,', error_count, 'errors.',
                         'Time (s):', duration_s)


if __name__ == "__main__":
    main()
