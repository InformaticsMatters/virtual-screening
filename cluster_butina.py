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

import argparse, time

from rdkit import DataStructs, rdBase
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.ML.Cluster import Butina

import utils
import rdkit_utils
from dm_job_utilities.dm_log import DmLog

descriptors = {
    'maccs':   lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan2': lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024),
    'morgan3': lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, 1024),
    'rdkit':   lambda m: AllChem.RDKFingerprint(m)
}

metrics = {
    'braunblanquet': DataStructs.BulkBraunBlanquetSimilarity,
    'cosine':        DataStructs.BulkCosineSimilarity,
    'dice':          DataStructs.BulkDiceSimilarity,
    'kulczynski':    DataStructs.BulkKulczynskiSimilarity,
    'mcconnaughey':  DataStructs.BulkMcConnaugheySimilarity,
    'rogotgoldberg': DataStructs.BulkRogotGoldbergSimilarity,
    'russel':        DataStructs.BulkRusselSimilarity,
    'sokal':         DataStructs.BulkSokalSimilarity,
    'tanimoto':      DataStructs.BulkTanimotoSimilarity
}

# start field name definitions #########################################

field_Cluster = "Cluster"

# functions ############################################################


def cluster_fps(fps, metric, cutoff):

    # first generate the distance matrix:
    dists = []
    # dist is the part of the distance matrix below the diagonal as an array:
    # 1.0, 2.0, 2.1, 3.0, 3.1, 3.2 ...
    nfps = len(fps)
    matrix = []
    for i in range(1, nfps):
        sims = metric(fps[i], fps[:i])
        dists.extend([1-x for x in sims])
        matrix.append(sims)

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs, dists, matrix


def clusters_to_map(clusters):
    d = {}
    i = 0
    for c in clusters:
        for id in c:
            d[id] = i
        i += 1
    return d


def fetch_score(idx, mols, field, descending):
    if descending:
        return 0 - mols[idx].GetDoubleProp(field)
    else:
        return mols[idx].GetDoubleProp(field)


def select_diverse_subset(mols, clusters, distances, count, field, descending, score):
    t0 = time.time()
    total = len(mols)
    num_clusters = len(clusters)
    picked_list = []
    clusters_list = []
    for i in range(0, num_clusters):
        picked_list.append([])
        if field:
            filtered_by_value = [x for x in clusters[i] if mols[x].HasProp(field)]
            sorted_by_value = sorted(filtered_by_value, key=lambda idx: fetch_score(idx, mols, field, descending))
            clusters_list.append(sorted_by_value)
        else:
            all_records = [x for x in clusters[i]]
            clusters_list.append(all_records)

    total_iter = 0
    cluster_iter = 0
    picked_count = 0

    while total_iter < total and picked_count < count:
        cluster_num = total_iter % num_clusters
        clus = clusters_list[cluster_num]
        pick = picked_list[cluster_num]
        if len(clus) > 0:
            # remove that item from the cluster so that it's not tried again
            mol_index = clus.pop(0)
            if len(pick) == 0:  # first time for this cluster
                pick.append(mol_index)
                picked_count += 1
                cluster_iter += 1
                # utils.log("Cluster", cluster_num, "initialised with molecule", mol_index)
            else:
                closest_dist = get_closest_distance(distances, mol_index, pick)
                if not score or closest_dist < score:
                    pick.append(mol_index)
                    picked_count += 1
                    cluster_iter += 1
                    # utils.log("Cluster", cluster_num, "added", mol_index, "with score", closestDist)
                # else:
                    # utils.log("Cluster", cluster_num, "discarded", mol_index, "with score", closestDist)
        else:  # cluster has been exhausted
            cluster_iter += 1

        total_iter += 1
    t1 = time.time()
    DmLog.emit_event("Picked {} molecules using {} iterations in {:.1f}s".format(picked_count, total_iter, t1 - t0))
    return picked_list


def get_distance(idx1, idx2, distances):
    idx = 0
    for i in range(1, idx1):
        idx += i
    idx += idx2
    d = distances[idx]
    return d


def get_closest_distance(distances, mol_idx, compare_to):
    best = 0
    for i in compare_to:
        d = get_distance(mol_idx, i, distances)
        if best < d:
            best = d
    return best


def execute(input, output, descriptor, metric, threshold, fragment_method, output_fragment, num, field, descending, exclude,
    delimiter=None, id_column=None, mol_column=0, omit_fields=False,
    read_header=False, write_header=False, read_records=50):

    # create reader
    calc_prop_names = [field_Cluster]
    reader = rdkit_utils.create_reader(input, id_column=id_column, mol_column=mol_column, read_records=read_records,
                                       read_header=read_header, delimiter=delimiter)
    extra_field_names = reader.get_extra_field_names()

    # create writer
    utils.expand_path(output)
    writer = rdkit_utils.create_writer(output,
                                       extra_field_names=extra_field_names,
                                       calc_prop_names=calc_prop_names,
                                       delimiter=delimiter,
                                       id_column=id_column, mol_column=mol_column)

    id_col_type, id_col_value = utils.is_type(id_column, int)

    # fragment and generate fingerprints
    mols = []  # the RDKit molecules
    data = []  # contains tuples of (id, smiles, props_dict) for each molecule
    fps = []   # the fingerprints for each molecule
    t0 = time.time()
    num_errs = rdkit_utils.fragmentAndFingerprint(reader, mols, data, fps, descriptor,
                                                  fragmentMethod=fragment_method, outputFragment=output_fragment)
    t1 = time.time()
    DmLog.emit_event("Read and fingerprinted {} molecules in {:.1f}s".format(len(mols), t1 - t0))
    if num_errs:
        DmLog.emit_event("Encountered {} errors fingerprinting molecules".format(num_errs))

    # do clustering
    t0 = time.time()
    clusters, dists, matrix = cluster_fps(fps, metric, 1.0 - threshold)
    t1 = time.time()
    DmLog.emit_event("Found {} clusters in {:.1f}s".format(len(clusters), t1 - t0))

    # generate diverse subset if specified
    # Note: max_min_picker.py is a much more scalable alternative
    if num:
        final_clusters = select_diverse_subset(mols, clusters, dists, num, field, descending, exclude)
    else:
        final_clusters = clusters

    # write the results
    lookup = clusters_to_map(final_clusters)
    i = 0
    result_count = 0
    for mol in mols:

        if result_count == 0 and write_header:
            headers = rdkit_utils.generate_headers(
                id_col_type,
                id_col_value,
                reader.get_mol_field_name(),
                reader.field_names,
                calc_prop_names,
                omit_fields)

            writer.write_header(headers)

        if i in lookup:
            cluster = lookup[i]
            writer.write(data[i][1], mol, data[i][0], data[i][2], (cluster,))
            result_count += 1
        i += 1
    DmLog.emit_event("Output {} molecules".format(result_count))
    DmLog.emit_cost(result_count)

    return len(clusters)


# start main execution ######################################################

def main():

    # Examples:
    #   python -m cluster_butina -i data/100.smi -o clustered.smi --id-column 1 -d tab --write-header -t 0.3

    # command line args definitions #########################################
    parser = argparse.ArgumentParser(description='RDKit Butina Cluster')
    parser.add_argument('-i', '--input', required=True, help="File with molecules to cluster (.sdf or .smi)")
    parser.add_argument('-o', '--output', required=True, help="Output file (.sdf or .smi)")

    parser.add_argument('-k', '--omit-fields', action='store_true',
                        help="Don't include fields from the input in the output")

    # to pass tab as the delimiter specify it as $'\t' or use one of the symbolic names 'comma', 'tab', 'space' or 'pipe'
    parser.add_argument('-d', '--delimiter', help="Delimiter when using SMILES")
    parser.add_argument('--id-column', help="Column for name field (zero based integer for .smi, text for SDF)")
    parser.add_argument('--mol-column', type=int, default=0,
                        help="Column index for molecule when using delineated text formats (zero based integer)")
    parser.add_argument('--read-header', action='store_true',
                        help="Read a header line with the field names when reading .smi or .txt")
    parser.add_argument('--write-header', action='store_true', help='Write a header line when writing .smi or .txt')
    parser.add_argument('--read-records', default=100, type=int,
                        help="Read this many records to determine the fields that are present")

    parser.add_argument('-t', '--threshold', type=float, default=0.7,
                        help='similarity clustering threshold (1.0 means identical)')
    parser.add_argument('--descriptor', type=str.lower, choices=list(descriptors.keys()), default='rdkit',
                        help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric', type=str.lower, choices=list(metrics.keys()), default='tanimoto',
                        help='similarity metric (default tanimoto)')
    parser.add_argument('-n', '--num', type=int, help='maximum number to pick for diverse subset selection')
    parser.add_argument('-e', '--exclude', type=float,
                        help='threshold for excluding structures in diverse subset selection (1.0 means identical)')
    parser.add_argument('--fragment-method', choices=['hac', 'mw'], default='hac',
                        help='How to find biggest fragment (hac: biggest by heavy atom count, mw: biggest by mol weight)')
    parser.add_argument('--output-fragment', action='store_true',
                        help='Output the biggest fragment rather than the original molecule')
    parser.add_argument('-f', '--field', help='field to use to optimise diverse subset selection')
    parser.add_argument('-a', '--ascending', action='store_true', help='Pick lowest value specified by the --field option')

    args = parser.parse_args()
    DmLog.emit_event("Cluster Butina Args: ", args)

    delimiter = utils.read_delimiter(args.delimiter)
    descriptor = descriptors[args.descriptor]
    metric = metrics[args.metric]

    if args.field and not args.num:
        raise ValueError('--num argument must be specified for diverse subset selection')

    num_clusters = execute(args.input, args.output, descriptor, metric, args.threshold, args.fragment_method,
                           args.output_fragment, args.num, args.field, not args.ascending, args.exclude,
                           omit_fields=args.omit_fields, delimiter=delimiter, id_column=args.id_column, mol_column=args.mol_column,
                           read_header=args.read_header, write_header=args.write_header,
                           read_records=args.read_records)


if __name__ == "__main__":
    main()

