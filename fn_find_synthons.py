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

# Cypher query that is executed:
# MATCH (fa:F2 {smiles: $smiles})-[e:FRAG*]->(f:F2) RETURN e
#


"""
Use the fragment network to find synthons of particular molecule
"""

import os, argparse, time
from rdkit import RDLogger
from neo4j import GraphDatabase
import utils
from standardize_molecule import standardize_to_noniso_smiles


RDLogger.logger().setLevel(RDLogger.ERROR)


def process(smiles, outfile, neo4j_server, neo4j_username, neo4j_password, standardize=True, report_hits=False):

    if standardize:
        std_smi, mol = standardize_to_noniso_smiles(smiles)
        utils.log('Standardized molecule:', smiles, ' -> ', std_smi)
        smiles = std_smi

    driver = GraphDatabase.driver(neo4j_server, auth=(neo4j_username, neo4j_password))
    with driver.session() as session:
        results = session.read_transaction(run_query, smiles)

    if report_hits:
        for result in results:
            utils.log(result)

    if outfile:
        with open(outfile, 'wt') as out:
            for result in results:
                out.write(result + '\n')

    return len(results)


def run_query(tx, smiles):

    query = "MATCH (fa:F2 {smiles: $smiles})-[e:FRAG*]->(f:F2) RETURN e"

    t0 = time.time()
    result = tx.run(query, smiles=smiles)
    t1 = time.time()
    synthons = set()
    for record in result:
        edges = record['e']
        for edge in edges:
            s = edge['label']
            tokens = s.split('|')
            synthons.add(tokens[1])
            synthons.add(tokens[4])

    utils.log_dm_event("Search Found", len(synthons), "synthons in", (t1 - t0), 'secs')
    return synthons


def main():
    # Example usage:
    # ./fn_find_synthons.py -s 'c1coc(CNc2nc3ccccc3[nH]2)c1' --server "bolt://localhost:7687" \
    #    --username <username> --password <password> --report-hits

    parser = argparse.ArgumentParser(description='Fragnet find synthons')
    parser.add_argument('-s', '--smiles', required=True, help='Query SMILES')
    parser.add_argument('-o', '--outfile', help='Output file')
    parser.add_argument('--no-standardize', action='store_true', help='The molecules do not need standardizing')
    parser.add_argument('--report-hits', action='store_true', help='Write the results to stdout')
    parser.add_argument('--server', help='Neo4j database server')
    parser.add_argument('--username', help='Neo4j username')
    parser.add_argument('--password', help='Neo4j password')

    args = parser.parse_args()
    utils.log_dm_event("FindSynthons Args: ", args)

    if args.server:
        neo4j_server = args.server
    else:
        neo4j_server = os.getenv('NEO4J_SERVER')
    if not neo4j_server:
        utils.log('WARNING: no database server URL provided')

    if args.username:
        neo4j_username = args.username
    else:
        neo4j_username = os.getenv('NEO4J_USERNAME')
    if not neo4j_username:
        utils.log('WARNING: no database username provided')

    if args.password:
        neo4j_password = args.password
    else:
        neo4j_password = os.getenv('NEO4J_PASSWORD')
    if not neo4j_password:
        utils.log('WARNING: no database password provided')

    # this does the processing
    count = process(args.smiles, args.outfile,
                    neo4j_server, neo4j_username, neo4j_password,
                    standardize=not args.no_standardize, report_hits=args.report_hits)

    utils.log_dm_event('Search found {} synthons.'.format(count))
    if count:
        utils.log_dm_cost(1)


if __name__ == "__main__":
    main()
